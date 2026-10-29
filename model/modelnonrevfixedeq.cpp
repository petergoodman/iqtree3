/*
 * modelnonrevfixedeq.cpp
 *
 *  Created on: 20/10/2026
 *      Author: Peter Goodman
 */

#include "modelnonrevfixedeq.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>

#include <iqtree_config.h>

#include "nclextra/modelsblock.h"
#include "partitionmodel.h"
#include "tree/iqtree.h"
#include "tree/phylotree.h"
#include "utils/eigendecomposition.h"

using namespace Eigen;

const char *const ModelNonrevFixedEq::NAME = "NQC";

static const char *const CHART = "jump-chain log-ratio";
static const int CHART_VERSION = 1;
static const char *const DERIVATIVE = "forward difference, step 1e-4 max(1 |x|), backward at the upper bound";
static const char *const STATES = "ARNDCQEGHILKMFPSTWYV";

static string number17(double v) {
    ostringstream ss;
    ss.precision(17);
    ss << v;
    return ss.str();
}

static string numbers17(const double *v, int n, const char *separator = ",") {
    string s;
    for (int i = 0; i < n; i++)
        s += (i ? separator : "") + number17(v[i]);
    return s;
}

static string integers(const std::vector<int> &v) {
    string s;
    for (size_t i = 0; i < v.size(); i++)
        s += (i ? "," : "") + convertIntToString(v[i]);
    return s;
}

/** a score, or "unknown" for NaN, which marks a score an older checkpoint did not hold */
static string score17(double v) {
    // saveCheckpoint writes only finite scores
    return std::isfinite(v) ? number17(v) : "unknown";
}

/** text on one line, each run of whitespace one space: a target defined over several lines of an
    --mdef file keeps its line breaks */
static string oneLine(const string &text) {
    string s;
    for (char c : text) {
        bool space = isspace((unsigned char)c);
        if (!space)
            s += c;
        else if (!s.empty() && s.back() != ' ')
            s += ' ';
    }
    if (!s.empty() && s.back() == ' ')
        s.pop_back();
    return s;
}

bool ModelNonrevFixedEq::isName(const string &model_name) {
    string upper_name = model_name;
    transform(upper_name.begin(), upper_name.end(), upper_name.begin(), ::toupper);
    return upper_name == NAME;
}

void ModelNonrevFixedEq::refuseAsMixtureComponent(ModelSubst *model) {
    if (dynamic_cast<ModelNonrevFixedEq*>(model))
        outError("NQC cannot be a component of a mixture model (MIX{...} or +FMIX)");
}

ModelNonrevFixedEq::ModelNonrevFixedEq(string model_params, StateFreqType freq, string freq_params,
                                       PhyloTree *tree, ModelsBlock *models_block)
    : ModelProtein(tree, models_block)
{
    ASSERT(num_states == 20);
    constructed = false;
    models = models_block;
    build_failed = false;
    fits = accepted_fit = 0;
    evaluations = 0;
    fit_starting = false;
    gradient_printed = false;
    start_lnl = optimizer_lnl = rescored_lnl = std::numeric_limits<double>::quiet_NaN();
    at_lower = at_upper = outside = 0;
    residuals = fixedeq::Diagnostics();
    refuseUnsupported(model_params, freq, freq_params);

    Params &params = Params::getInstance();
    string error;
    if (fixedeq::parseTarget(freq_params, num_states, max(params.min_state_freq, ZERO_FREQ), target,
                             target_raw_sum, error) != fixedeq::Status::ok)
        outError("NQC: invalid target +F{" + freq_params + "}: " + error);
    target_text = freq_params;
    lower = log(params.nqc_domain_lo);
    upper = log(params.nqc_domain_hi);

    MatrixXd seed;
    readSeed(seed);
    VectorXd z;
    if (fixedeq::inverse(seed, ref, z) != fixedeq::Status::ok)
        outError("NQC: the starting matrix has a rate that is zero, negative or not finite");
    trial = accepted = z;
    Q.resize(num_states, num_states);
    nu.resize(num_states);
    freq_type = FREQ_USER_DEFINED;
    num_params = (int)z.size();
    name = NAME;
    full_name = "Non-reversible model with a fixed stationary distribution (NQC)";
    constructed = true;
    decomposeRateMatrix();
    if (build_failed)
        outError("NQC: the starting coordinates do not give a valid rate matrix");
    measureState();
    if (verbose_mode >= VB_MED) {
        cout << "NQC seed Q after rescaling (17 significant digits):" << endl;
        printMatrix(cout, rate_matrix);
    }
}

void ModelNonrevFixedEq::refuseUnsupported(const string &model_params, StateFreqType freq,
                                           const string &freq_params) {
    Params &params = Params::getInstance();
    const char *why = nullptr;
    // runModelFinder sets dating_mf, although named for dating, before it builds any candidate,
    // so it is the run-time sign that ModelFinder or MixtureFinder is testing NQC. It and the
    // site-specific frequency check come first, because a candidate's string, and the string
    // the site-specific path builds each model from, lack the target.
    if (params.dating_mf)
        why = "cannot be tested by ModelFinder; give it with -m";
    else if (params.site_freq_file || params.tree_freq_file)
        why = "is not supported with site-specific frequency models";
    else if (!model_params.empty())
        why = "takes no parameters in braces after its name; give its target as +F{...}";
    else if (freq != FREQ_USER_DEFINED || freq_params.empty())
        why = "needs its target stationary distribution, as +F{p1,...,p20} or as +F<name> defined in an "
              "--mdef file; for frequencies counted (+F) or estimated (+FO) from the data use NONREV";
    else if (params.optimize_from_given_params)
        why = "is not supported with -optfromgiven";
    else if (params.optimize_model_rate_joint)
        why = "is not supported with -jointopt";
    else if (params.model_name.find("+T") != string::npos || params.optimize_params_use_hmm)
        why = "is not supported in tree mixtures or HMM models";
    // a partition's own model string, as a charset gives it, never reaches params.model_name
    else if (params.model_name.find("+H") != string::npos ||
             (phylo_tree && phylo_tree->aln && phylo_tree->aln->model_name.find("+H") != string::npos))
        why = "is not supported with heterotachy (+H) models";
    else if (params.optimize_linked_gtr)
        why = "is not supported with --link-exchange-rates";
    else if (params.matrix_exp_technique == MET_EIGEN_DECOMPOSITION)
        why = "is not supported with --eigen";
    else if (params.alisim_active)
        why = "cannot simulate alignments; simulate from its fitted matrix instead";
    // after the AliSim check, because AliSim also sets partition_file
    else if (params.partition_file && params.model_joint.empty() && !params.link_model)
        why = "is not supported with unlinked partitions (-p, -q or -S without --model-joint or --link-model)";
    // cached transition matrices, keyed by branch length, would outlive a change of Q
    else if (params.store_trans_matrix)
        why = "is not supported with -mstore";
    // --runs restores the best run into a tree that cannot compute a likelihood, which the
    // export's diagnostic needs
    else if (params.num_runs > 1)
        why = "is not supported with --runs; for independent runs, run IQ-TREE once per -seed";
    if (why)
        outError(string("NQC ") + why);
}

void ModelNonrevFixedEq::readSeed(MatrixXd &seed) {
    // Read with tree adaptation off, as ModelProtein::init's NONREV branch does, so that a rooted
    // input tree keeps its root (decision 031). No rescaling: both seed routes are scale-invariant.
    NxsModel *lg = models->findModel("LG");
    ASSERT(lg);
    readParametersString(lg->description, false);
    MatrixXd R = MatrixXd::Zero(num_states, num_states);
    for (int i = 0, k = 0; i < num_states; i++)
        for (int j = i + 1; j < num_states; j++, k++)
            R(i, j) = R(j, i) = rates[k];
    if (fixedeq::reversibleSeed(R, target, seed) != fixedeq::Status::ok)
        outError("NQC: LG cannot be rebuilt at the target");
    ref = fixedeq::references(seed);
    start_text = "LG rebuilt at the target";

    const char *init = Params::getInstance().model_name_init;
    if (init) {
        NxsModel *nxs = models->findModel(init);
        if (nxs)
            readParametersString(nxs->description, false);
        else
            readParameters(init, false);
        start_text = string("--init-model ") + (nxs ? "model " : "file ") + init +
                     (is_reversible ? ", rebuilt at the target" : ", its jump chain transferred to the target");
        fixedeq::Status status;
        if (is_reversible) {
            for (int i = 0, k = 0; i < num_states; i++)
                for (int j = i + 1; j < num_states; j++, k++)
                    R(i, j) = R(j, i) = rates[k];
            status = fixedeq::reversibleSeed(R, target, seed);
        } else {
            MatrixXd Q0(num_states, num_states);
            for (int i = 0, k = 0; i < num_states; i++) {
                for (int j = 0; j < num_states; j++)
                    if (j != i)
                        Q0(i, j) = rates[k++];
                Q0(i, i) = 0.0;
                Q0(i, i) = -Q0.row(i).sum();
            }
            status = fixedeq::transferSeed(Q0, target, seed);
        }
        if (status != fixedeq::Status::ok)
            outError(string("NQC: the --init-model start cannot be moved to the target: ") + fixedeq::statusText(status));
    }

    memcpy(state_freq, target.data(), sizeof(double) * num_states);
    if (is_reversible) {
        // also roots an unrooted tree and keeps a rooted one
        setReversible(false);
    } else if (phylo_tree && !phylo_tree->rooted && phylo_tree->root) {
        // the rooting setReversible(false) would do, which a model read as non-reversible skips
        if (verbose_mode >= VB_MED)
            cout << "Converting unrooted to rooted tree..." << endl;
        phylo_tree->convertToRooted();
    }
}

void ModelNonrevFixedEq::setVariables(double *variables) {
    int ndim = getNDim();
    if (ndim > 0)
        memcpy(variables + 1, accepted.data(), ndim * sizeof(double));
}

bool ModelNonrevFixedEq::getVariables(double *variables) {
    int ndim = getNDim();
    bool changed = false;
    for (int i = 0; i < ndim; i++)
        changed |= (trial[i] != variables[i + 1]);
    if (ndim > 0)
        memcpy(trial.data(), variables + 1, ndim * sizeof(double));
    return changed;
}

void ModelNonrevFixedEq::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    int ndim = getNDim();
    for (int i = 1; i <= ndim; i++) {
        lower_bound[i] = lower;
        upper_bound[i] = upper;
        bound_check[i] = false;
    }
}

double ModelNonrevFixedEq::targetFunk(double x[]) {
    evaluations++;
    bool changed = getVariables(x);
    if (changed) {
        decomposeRateMatrix();
        ASSERT(phylo_tree);
        if (!build_failed)
            phylo_tree->clearAllPartialLH();
    }
    // The base routine's min_state_freq guard is not repeated: pi* was checked against it.
    if (build_failed)
        return 1.0e+30;
    return -phylo_tree->computeLikelihood();
}

double ModelNonrevFixedEq::derivativeFunk(double x[], double dfx[]) {
    return scaledDerivative(*this, x, dfx);
}

double ModelNonrevFixedEq::linkedDerivative(Optimization &objective, double x[], double dfx[]) {
    return scaledDerivative(objective, x, dfx);
}

double ModelNonrevFixedEq::scaledDerivative(Optimization &objective, double x[], double dfx[]) {
    // the objective is this model in a single-model fit and the PartitionModel in a linked one
    int ndim = getNDim();
    double *h = new double[ndim + 1];
    double fx = objective.targetFunk(x);
    for (int dim = 1; dim <= ndim; dim++) {
        double temp = x[dim];
        x[dim] = temp + fixedeq::scaledStep(temp, upper);
        h[dim] = x[dim] - temp;
        dfx[dim] = objective.targetFunk(x);
        x[dim] = temp;
    }
    for (int dim = 1; dim <= ndim; dim++)
        dfx[dim] = (dfx[dim] - fx) / h[dim];
    if (fit_starting) {
        fit_starting = false;
        start_lnl = -fx;
        if (verbose_mode >= VB_MAX && !gradient_printed) {
            gradient_printed = true;
            cout << "NQC first gradient at -lnL " << number17(fx)
                 << " (coordinate index, coordinate, step, derivative of -lnL):" << endl;
            for (int dim = 1; dim <= ndim; dim++)
                cout << "NQC gradient " << dim << " " << number17(x[dim]) << " " << number17(h[dim]) << " "
                     << number17(dfx[dim]) << endl;
        }
    }
    delete [] h;
    return fx;
}

double ModelNonrevFixedEq::optimizeParameters(double gradient_epsilon) {
    // ModelMarkov::optimizeParameters, statement for statement, without its unused second vector
    if (fixed_parameters)
        return 0.0;
    int ndim = getNDim();
    if (ndim == 0)
        return 0.0;
    if (verbose_mode >= VB_MAX)
        cout << "Optimizing " << name << " model parameters..." << endl;

    double *variables = new double[ndim + 1];
    double *upper_bound = new double[ndim + 1];
    double *lower_bound = new double[ndim + 1];
    bool *bound_check = new bool[ndim + 1];
    double score;

    for (int i = 0; i < num_states; i++)
        if (state_freq[i] > state_freq[highest_freq_state])
            highest_freq_state = i;

    setVariables(variables);
    setBounds(lower_bound, upper_bound, bound_check);
    fit_starting = true;
    score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(gradient_epsilon, TOL_RATE));
    double optimizer_score = score;
    bool changed = getVariables(variables);
    if (changed || score == -1.0e+30) {
        decomposeRateMatrix();
        phylo_tree->clearAllPartialLH();
        score = phylo_tree->computeLikelihood();
    }
    checkFit(optimizer_score);

    delete [] bound_check;
    delete [] lower_bound;
    delete [] upper_bound;
    delete [] variables;
    return score;
}

static const char *const BUILD_FAILED_AFTER_FIT =
    "NQC: the optimizer returned coordinates that do not give a valid rate matrix; the last valid coordinates are kept";

static bool sameBits(const double *a, const double *b, size_t n) {
    return memcmp(a, b, sizeof(double) * n) == 0;
}

void ModelNonrevFixedEq::checkFit(double optimizer_score) {
    // decision 009: re-score the returned point and report; nothing is restored. This runs only
    // outside parallel regions, because unlinked partitions are refused (decision 008).
    double lnl = rescore();
    if (build_failed)
        outWarning(BUILD_FAILED_AFTER_FIT);
    recordFit(start_lnl, optimizer_score, lnl);
    reportFit();
}

void ModelNonrevFixedEq::startLinkedFit(const std::vector<ModelSubst*> &linked) {
    // Checked before every linked fit, so the first check follows round 1's per-partition updates
    // and precedes any move of the shared matrix; those updates compute nothing wrong, since each
    // partition uses its own target (CODE_PLAN.md section 2.2). Chart, version and derivative
    // policy are constants of this class, so every NQC object shares them.
    double domain[2] = {lower, upper};
    for (ModelSubst *m : linked) {
        ModelNonrevFixedEq *other = dynamic_cast<ModelNonrevFixedEq*>(m);
        const char *why = nullptr;
        if (!other) {
            why = "a model that is not NQC";
        } else {
            double other_domain[2] = {other->lower, other->upper};
            if (!sameBits(target.data(), other->target.data(), num_states) ||
                !sameBits(state_freq, other->state_freq, num_states))
                why = "a different target";
            else if (other->ref != ref)
                why = "different reference destinations";
            else if (!sameBits(domain, other_domain, 2))
                why = "a different coordinate domain";
            else if (!sameBits(accepted.data(), other->accepted.data(), accepted.size()))
                why = "different coordinates";
        }
        if (why)
            outError(string("NQC: linked partitions must share one model, but partition ") +
                     (other ? other->phylo_tree->aln->name : string("?")) + " has " + why);
    }
    fit_starting = true;
}

void ModelNonrevFixedEq::endLinkedFit(const std::vector<ModelSubst*> &linked, double optimizer_score) {
    // decision 009 for a linked fit, outside targetFunk's parallel region. The re-score sums each
    // partition's lnL in the list's order, bitwise the sum PartitionModel::targetFunk forms, and
    // counts no evaluation; through targetFunk a failed build would score K x 1e30. Every object
    // records the fit, so all partition checkpoint structs agree, and this one, the representative,
    // reports it once, before the round's ASSERT.
    double lnl = 0.0;
    for (ModelSubst *m : linked)
        lnl += static_cast<ModelNonrevFixedEq*>(m)->rescore();
    if (build_failed)
        outWarning(BUILD_FAILED_AFTER_FIT);
    for (ModelSubst *m : linked)
        static_cast<ModelNonrevFixedEq*>(m)->recordFit(start_lnl, optimizer_score, lnl);
    reportFit();
}

double ModelNonrevFixedEq::rescore() {
    // Between fits IQ-TREE updates rates and branches without decomposing Q, so the decomposition's
    // guard on pi* is repeated here. A failed build keeps the last valid state, which is scored.
    if (memcmp(state_freq, target.data(), sizeof(double) * num_states) != 0)
        outError("NQC: the stationary distribution was changed between fits");
    decomposeRateMatrix();
    phylo_tree->clearAllPartialLH();
    return phylo_tree->computeLikelihood();
}

void ModelNonrevFixedEq::recordFit(double start, double optimizer_score, double rescored) {
    fits++;
    accepted_fit = fits;
    start_lnl = start;
    optimizer_lnl = optimizer_score;
    rescored_lnl = rescored;
    measureState();
}

void ModelNonrevFixedEq::reportFit() {
    if (verbose_mode >= VB_MED) {
        cout << "NQC fit " << fits << ": start lnL " << number17(start_lnl) << ", optimizer lnL "
             << number17(optimizer_lnl) << ", re-scored lnL " << number17(rescored_lnl) << endl;
        writeState(cout);
    }
    if (rescored_lnl < start_lnl)
        outWarning("NQC fit " + convertIntToString(fits) + ": the re-scored lnL is below the start by " +
                   number17(start_lnl - rescored_lnl) + " (reported, not restored)");
}

void ModelNonrevFixedEq::measureState() {
    // derived from the accepted coordinates and Q, so recomputed rather than restored
    at_lower = at_upper = outside = 0;
    for (int i = 0; i < accepted.size(); i++) {
        at_lower += (accepted[i] == lower);
        at_upper += (accepted[i] == upper);
        outside += (accepted[i] < lower || accepted[i] > upper);
    }
    Map<const Matrix<double, Dynamic, Dynamic, RowMajor> > q(rate_matrix, num_states, num_states);
    residuals = fixedeq::diagnostics(q, target, nu);
}

void ModelNonrevFixedEq::writeState(ostream &out) {
    out << "NQC coordinates at the lower bound " << at_lower << ", at the upper bound " << at_upper
        << ", outside the domain " << outside << endl;
    out << "NQC residuals: |pi*Q| " << number17(residuals.stationarity) << ", row sums "
        << number17(residuals.row_sums) << ", |mean rate - 1| " << number17(residuals.mean_rate_error)
        << ", nu " << number17(residuals.nu_residual) << ", minimum rate " << number17(residuals.min_rate) << endl;
}

void ModelNonrevFixedEq::decomposeRateMatrixNonrev() {
    if (!constructed) {
        ModelMarkov::decomposeRateMatrixNonrev();
        return;
    }
    fixedeq::Status status = fixedeq::build(trial, target, ref, Q, nu);
    if (status != fixedeq::Status::ok) {
        // rates[], the eigensystem and the accepted coordinates stay as they were
        build_failed = true;
        if (verbose_mode >= VB_DEBUG) {
            // in a linked fit this runs inside PartitionModel::targetFunk's parallel loop
#ifdef _OPENMP
#pragma omp critical
#endif
            cout << "NQC: proposal rejected: " << fixedeq::statusText(status) << endl;
        }
        return;
    }
    build_failed = false;
    accepted = trial;
    for (int i = 0, k = 0; i < num_states; i++)
        for (int j = 0; j < num_states; j++)
            if (j != i)
                rates[k++] = Q(i, j);
    ModelMarkov::decomposeRateMatrixNonrev();
    if (memcmp(state_freq, target.data(), sizeof(double) * num_states) != 0)
        outError("NQC: the stationary distribution was changed during decomposition");
}

void ModelNonrevFixedEq::getStateFrequency(double *freq, int mixture) {
    if (!constructed) {
        ModelMarkov::getStateFrequency(freq, mixture);
        return;
    }
    memcpy(freq, state_freq, sizeof(double) * num_states);
}

void ModelNonrevFixedEq::setStateFrequency(double *freq) {
    if (constructed && memcmp(freq, target.data(), sizeof(double) * num_states) != 0)
        outError("NQC: the stationary distribution is fixed to its target and cannot be set");
    ModelMarkov::setStateFrequency(freq);
}

void ModelNonrevFixedEq::adaptStateFrequency(double *freq) {
    if (!constructed) {
        ModelMarkov::adaptStateFrequency(freq);
        return;
    }
    if (memcmp(freq, target.data(), sizeof(double) * num_states) != 0)
        outError("NQC: the stationary distribution is fixed to its target and cannot be changed");
}

void ModelNonrevFixedEq::setRateMatrix(double *rate_mat) {
    if (constructed)
        outError("NQC: the rate matrix is set only through the model's own coordinates");
    ModelMarkov::setRateMatrix(rate_mat);
}

void ModelNonrevFixedEq::setFullRateMatrix(double *rate_mat, double *freq) {
    if (constructed)
        outError("NQC: the rate matrix is set only through the model's own coordinates");
    ModelMarkov::setFullRateMatrix(rate_mat, freq);
}

void ModelNonrevFixedEq::startCheckpoint() {
    checkpoint->startStruct("ModelNonrevFixedEq");
}

void ModelNonrevFixedEq::saveCheckpoint() {
    // 17-digit strings rather than CKP_SAVE's 10 digits (decision 003)
    double domain[2] = {lower, upper};
    // the seven diagnostics in fixedeq::Diagnostics order
    double diagnostics[7] = {residuals.stationarity, residuals.row_sums, residuals.mean_rate_error,
                             residuals.min_rate, residuals.min_flux, residuals.nu_residual,
                             residuals.exit_rate_mismatch};
    double scores[3] = {start_lnl, optimizer_lnl, rescored_lnl};
    const char *score_keys[3] = {"start_score", "optimizer_score", "accepted_score"};
    startCheckpoint();
    checkpoint->put("chart", string(CHART));
    checkpoint->put("chart_version", CHART_VERSION);
    checkpoint->put("target", numbers17(state_freq, num_states));
    checkpoint->put("references", integers(ref));
    checkpoint->put("domain", numbers17(domain, 2));
    checkpoint->put("derivative", string(DERIVATIVE));
    checkpoint->put("coordinates", numbers17(accepted.data(), (int)accepted.size()));
    checkpoint->putBool("fixed", fixed_parameters);
    checkpoint->put("fits", fits);
    checkpoint->put("evaluations", evaluations);
    checkpoint->put("accepted_fit", accepted_fit);
    for (int i = 0; i < 3; i++)
        if (std::isfinite(scores[i]))
            checkpoint->put(score_keys[i], number17(scores[i]));
    checkpoint->put("bounds", integers({at_lower, at_upper, outside}));
    checkpoint->put("residuals", numbers17(diagnostics, 7));
    endCheckpoint();
    ModelProtein::saveCheckpoint();
}

void ModelNonrevFixedEq::restoreCheckpoint() {
    // The base restore decomposes from the coordinates held now; this restore then rebuilds from
    // the saved ones. Absent keys leave the model as it is, as IQ-TREE's own restores do.
    ModelProtein::restoreCheckpoint();
    double domain[2] = {lower, upper};
    double *scores[3] = {&start_lnl, &optimizer_lnl, &rescored_lnl};
    const char *score_keys[3] = {"start_score", "optimizer_score", "accepted_score"};
    string coordinates, chart, saved_target, saved_ref, saved_domain;
    int chart_version = 0, saved_fits = 0;
    int64_t saved_evaluations = 0;
    startCheckpoint();
    bool found = checkpoint->getString("coordinates", coordinates);
    if (found) {
        checkpoint->getString("chart", chart);
        checkpoint->get("chart_version", chart_version);
        checkpoint->getString("target", saved_target);
        checkpoint->getString("references", saved_ref);
        checkpoint->getString("domain", saved_domain);
        if (chart != CHART || chart_version != CHART_VERSION || saved_target != numbers17(target.data(), num_states) ||
            saved_ref != integers(ref) || saved_domain != numbers17(domain, 2))
            outError("NQC: the checkpoint was written with a different target, references, chart or domain");
        // counters keep the larger count, so that neither the +I+G path's restore of its best start
        // nor a resumed run winds them back
        checkpoint->get("fits", saved_fits);
        checkpoint->get("evaluations", saved_evaluations);
        fits = max(fits, saved_fits);
        evaluations = max(evaluations, saved_evaluations);
        accepted_fit = 0;
        checkpoint->get("accepted_fit", accepted_fit);
        for (int i = 0; i < 3; i++) {
            string score;
            *scores[i] = checkpoint->getString(score_keys[i], score) ? convert_double(score.c_str())
                                                                    : std::numeric_limits<double>::quiet_NaN();
        }
    }
    endCheckpoint();
    if (!found)
        return;
    DoubleVector z;
    convert_double_vec(coordinates.c_str(), z);
    ASSERT(z.size() == (size_t)trial.size());
    memcpy(trial.data(), z.data(), sizeof(double) * z.size());
    decomposeRateMatrix();
    if (build_failed)
        outError("NQC: the checkpointed coordinates do not give a valid rate matrix");
    measureState();
    if (phylo_tree)
        phylo_tree->clearAllPartialLH();
}

void ModelNonrevFixedEq::printMatrix(ostream &out, const double *m) {
    for (int i = 0; i < num_states; i++)
        out << numbers17(m + i * num_states, num_states) << endl;
}

void ModelNonrevFixedEq::reportRun(IQTree &tree, ostream &out) {
    // Partitioned NQC is always linked, since unlinked partitions are refused (decision 008), so a
    // partitioned run gets one export and one block per linked model, from its representative.
    if (!tree.isSuperTree()) {
        ModelNonrevFixedEq *model = dynamic_cast<ModelNonrevFixedEq*>(tree.getModel());
        if (model)
            model->exportRun(tree, *model, {model}, 0, out);
        return;
    }
    PartitionModel *partition_model = dynamic_cast<PartitionModel*>(tree.getModelFactory());
    if (!partition_model)
        return;
    PhyloSuperTree &super_tree = (PhyloSuperTree&)tree;
    for (auto &linked : partition_model->linked_models) {
        ModelNonrevFixedEq *model = dynamic_cast<ModelNonrevFixedEq*>(linked.second);
        if (!model)
            continue;
        // built only once an NQC model is found, since at -v building the order prints it, and a
        // run without NQC must print nothing new (PLAN.md constraint 4)
        if (super_tree.part_order.empty())
            super_tree.computePartitionOrder();
        std::vector<ModelNonrevFixedEq*> objects;
        for (int j = 0; j < (int)super_tree.size(); j++) {
            ModelSubst *part_model = super_tree.at(super_tree.part_order[j])->getModel();
            if (part_model->getName() == linked.first)
                objects.push_back(dynamic_cast<ModelNonrevFixedEq*>(part_model));
        }
        // PartitionModel::targetFunk evaluates the partitions whose model has the name of its model
        // member, a placeholder outside optimizeLinkedModels, which sets it as here. A linked model
        // stays fixed between its fits; it is reported unfixed, as PartitionModel saves and
        // reportPhyloAnalysis prints it, so that the export record matches its checkpoint.
        ModelSubst *saved_model = partition_model->model;
        partition_model->model = model;
        bool fixed = model->fixParameters(false);
        model->exportRun(tree, *partition_model, objects, (int)objects.size(), out);
        model->fixParameters(fixed);
        partition_model->model = saved_model;
    }
}

void ModelNonrevFixedEq::exportRun(PhyloTree &tree, Optimization &objective,
                                   const std::vector<ModelNonrevFixedEq*> &objects, int linked_partitions,
                                   ostream &out) {
    string prefix = Params::getInstance().out_prefix;
    string qmat_file = prefix + "." + NAME + ".qmat", info_file = prefix + "." + NAME + ".info";
    // the matrix first, so that nothing the diagnostic does can lose the fit
    writeQmat(qmat_file);
    out << endl << "NQC first-order diagnostic at export: " << 2 * accepted.size() + 2
        << " likelihood evaluations" << endl;
    ExportDiagnostic d = exportDiagnostic(tree, objective, objects);
    writeSummary(out, linked_partitions, d, qmat_file, info_file);
    writeExportInfo(info_file, linked_partitions, d, qmat_file);
}

void ModelNonrevFixedEq::writeQmat(const string &file) {
    // the matrix the likelihood used, normalized by the base decomposition, and pi*, in the -m FILE
    // format; no header, since the reader takes the first token's sign to choose a full matrix
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(file.c_str());
        for (int i = 0; i < num_states; i++)
            out << numbers17(rate_matrix + i * num_states, num_states, " ") << '\n';
        out << numbers17(state_freq, num_states, " ") << '\n';
        out.close();
    } catch (ios::failure &) {
        outError(ERR_WRITE_OUTPUT, file);
    }
}

ModelNonrevFixedEq::ExportDiagnostic ModelNonrevFixedEq::exportDiagnostic(
        PhyloTree &tree, Optimization &objective, const std::vector<ModelNonrevFixedEq*> &objects) {
    // Decision 009: central differences of -lnL through the objective the fits used, then the
    // accepted state rebuilt. Sized by the coordinates, since getNDim is 0 for a fixed model.
    const int n = (int)accepted.size();
    const VectorXd x0 = accepted;
    // the matrix the run ended with; every linked object holds the same one
    const std::vector<double> matrix(rate_matrix, rate_matrix + num_states * num_states);
    std::vector<int64_t> counts;
    for (ModelNonrevFixedEq *m : objects) {
        ASSERT(m);
        if (!sameBits(m->accepted.data(), x0.data(), n))
            outError("NQC: the linked partitions hold different coordinates at export");
        counts.push_back(m->evaluations);
        m->trial = x0;
        m->decomposeRateMatrix();
    }
    // every evaluation below starts from cleared partials; the -b report re-reads the tree without
    // reinitializing them, as optimizeConTree does after the same re-read
    tree.initializeAllPartialLh();
    tree.clearAllPartialLH();

    ExportDiagnostic d;
    d.g.resize(n);
    d.failed = 0;
    double *x = new double[n + 1];
    memcpy(x + 1, x0.data(), n * sizeof(double));
    d.f0 = objective.targetFunk(x);
    for (int k = 1; k <= n; k++) {
        double temp = x[k];
        double h = fixedeq::FIRST_ORDER_ETA * max(1.0, fabs(temp));
        x[k] = temp + h;
        double up_step = x[k] - temp;
        double up = objective.targetFunk(x);
        x[k] = temp - h;
        double down_step = temp - x[k];
        double down = objective.targetFunk(x);
        x[k] = temp;
        d.failed += (up >= 1.0e+30) + (down >= 1.0e+30);
        d.g[k - 1] = (up - down) / (up_step + down_step);
    }
    // the last probe moved a coordinate, so this rebuilds the accepted state
    double restored = objective.targetFunk(x);
    delete [] x;
    d.evaluations = 2 * n + 2;
    for (size_t i = 0; i < objects.size(); i++) {
        objects[i]->evaluations = counts[i];
        // the .iqtree prints this matrix, and -z evaluates with it
        ASSERT(memcmp(objects[i]->rate_matrix, matrix.data(), sizeof(double) * matrix.size()) == 0);
    }
    ASSERT(restored == d.f0);

    if (fixedeq::firstOrderTest(x0, d.g, lower, upper, d.f0, fixedeq::FIRST_ORDER_GTOL, d.test) !=
        fixedeq::Status::ok) {
        d.test = {std::numeric_limits<double>::quiet_NaN(), 0, 0, 0, 0, false};
        d.reason = "the differences are not finite";
    } else if (d.failed > 0) {
        d.reason = "a probe gave no valid matrix";
    } else if (d.test.outside > 0) {
        d.reason = "coordinates lie outside the domain";
    } else if (!d.test.passed) {
        d.reason = "the projected gradient fails the optimizer's test";
    } else {
        d.reason = "the projected gradient passes the optimizer's test";
    }
    d.label = (d.test.passed && d.failed == 0) ? "first-order stationary at the optimizer's tolerance"
                                                : "unclassified";
    return d;
}

void ModelNonrevFixedEq::writeSummary(ostream &out, int linked_partitions, const ExportDiagnostic &d,
                                      const string &qmat_file, const string &info_file) {
    string letters;
    for (int r : ref)
        letters += STATES[r];
    out << endl << "NQC summary of the final model (17 significant digits):" << endl;
    if (linked_partitions > 0)
        out << "NQC linked across " << linked_partitions << " partitions" << endl;
    out << "NQC target as given (sum " << number17(target_raw_sum) << "): " << target_text << endl;
    out << "NQC target normalized: " << numbers17(state_freq, num_states) << endl;
    out << "NQC reference destinations: " << letters << endl;
    out << "NQC coordinate domain: [" << number17(lower) << ", " << number17(upper) << "]" << endl;
    out << "NQC accepted fit " << (accepted_fit > 0 ? convertIntToString(accepted_fit) : string("unknown"))
        << " of " << fits << ": start lnL " << score17(start_lnl) << ", optimizer lnL " << score17(optimizer_lnl)
        << ", re-scored lnL " << score17(rescored_lnl) << endl;
    writeState(out);
    out << "NQC likelihood evaluations: " << evaluations << endl;
    out << "NQC stop: " << d.label << " (" << d.reason << "; largest scaled component "
        << number17(d.test.statistic) << " at coordinate " << d.test.coordinate + 1 << ", line "
        << number17(fixedeq::FIRST_ORDER_GTOL) << ")" << endl;
    out << "NQC fitted matrix written to " << qmat_file << "; re-use it with -m " << qmat_file << endl;
    out << "NQC export record written to " << info_file << endl;
}

void ModelNonrevFixedEq::writeExportInfo(const string &file, int linked_partitions, const ExportDiagnostic &d,
                                         const string &qmat_file) {
    // keys shared with the checkpoint are formatted as saveCheckpoint formats them
    double domain[2] = {lower, upper};
    double diagnostics[7] = {residuals.stationarity, residuals.row_sums, residuals.mean_rate_error,
                             residuals.min_rate, residuals.min_flux, residuals.nu_residual,
                             residuals.exit_rate_mismatch};
    string letters;
    for (int r : ref)
        letters += STATES[r];
    string command = "unknown";
    if (checkpoint) {
        checkpoint->startStruct("iqtree");
        checkpoint->getString("command", command);
        checkpoint->endStruct();
        // the record begins with the space before the first argument
        if (!command.empty() && command[0] == ' ')
            command.erase(0, 1);
    }
    ostringstream version;
    version << iqtree_VERSION_MAJOR << "." << iqtree_VERSION_MINOR << iqtree_VERSION_PATCH;
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(file.c_str());
        out << "# NQC export record of " << qmat_file << " (decisions 009 and 010)" << '\n'
            << "model: " << NAME << '\n'
            << "state_order: " << STATES << '\n'
            << "target_as_given: " << oneLine(target_text) << '\n'
            << "target_sum: " << number17(target_raw_sum) << '\n'
            << "target: " << numbers17(state_freq, num_states) << '\n'
            << "chart: " << CHART << '\n'
            << "chart_version: " << CHART_VERSION << '\n'
            << "references: " << integers(ref) << '\n'
            << "reference_letters: " << letters << '\n'
            << "domain: " << numbers17(domain, 2) << '\n'
            << "derivative: " << DERIVATIVE << '\n'
            << "start: " << start_text << '\n'
            << "random_seed: " << Params::getInstance().ran_seed << '\n'
            << "command: " << oneLine(command) << '\n'
            << "linked_partitions: " << linked_partitions << '\n'
            << "fixed: " << (fixed_parameters ? "true" : "false") << '\n'
            << "fits: " << fits << '\n'
            << "evaluations: " << evaluations << '\n'
            << "accepted_fit: " << accepted_fit << '\n'
            << "start_score: " << score17(start_lnl) << '\n'
            << "optimizer_score: " << score17(optimizer_lnl) << '\n'
            << "accepted_score: " << score17(rescored_lnl) << '\n'
            << "bounds: " << integers({at_lower, at_upper, outside}) << '\n'
            << "residuals: " << numbers17(diagnostics, 7) << '\n'
            << "residual_names: stationarity,row_sums,mean_rate_error,min_rate,min_flux,nu_residual,"
               "exit_rate_mismatch" << '\n'
            << "coordinates: " << numbers17(accepted.data(), (int)accepted.size()) << '\n'
            << "export_lnl: " << number17(-d.f0) << '\n'
            << "stop: " << d.label << '\n'
            << "stop_reason: " << d.reason << '\n'
            << "stop_statistic: " << number17(d.test.statistic) << '\n'
            << "stop_coordinate: " << d.test.coordinate + 1 << '\n'
            << "stop_line: " << number17(fixedeq::FIRST_ORDER_GTOL) << '\n'
            << "stop_step: central differences, step 1e-3 max(1 |x|), over the probes' displacements" << '\n'
            << "stop_failed_probes: " << d.failed << '\n'
            << "stop_gradient: " << numbers17(d.g.data(), (int)d.g.size()) << '\n'
            << "diagnostic_evaluations: " << d.evaluations << '\n'
            << "qmat: " << qmat_file << '\n'
            << "iqtree_version: " << version.str() << '\n';
        out.close();
    } catch (ios::failure &) {
        outError(ERR_WRITE_OUTPUT, file);
    }
}
