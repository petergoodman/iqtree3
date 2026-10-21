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
#include <sstream>

#include "nclextra/modelsblock.h"
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

static string numbers17(const double *v, int n) {
    string s;
    for (int i = 0; i < n; i++)
        s += (i ? "," : "") + number17(v[i]);
    return s;
}

static string integers(const std::vector<int> &v) {
    string s;
    for (size_t i = 0; i < v.size(); i++)
        s += (i ? "," : "") + convertIntToString(v[i]);
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
    fits = 0;
    fit_starting = false;
    gradient_printed = false;
    start_lnl = optimizer_lnl = rescored_lnl = 0.0;
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
    if (verbose_mode >= VB_MED) {
        cout << "NQC seed Q after rescaling (17 significant digits):" << endl;
        printMatrix(cout, rate_matrix);
    }
}

void ModelNonrevFixedEq::refuseUnsupported(const string &model_params, StateFreqType freq,
                                           const string &freq_params) {
    Params &params = Params::getInstance();
    const char *why = nullptr;
    if (!model_params.empty())
        why = "takes no parameters in braces after its name; give its target as +F{...}";
    else if (freq != FREQ_USER_DEFINED || freq_params.empty())
        why = "needs its target stationary distribution, as +F{p1,...,p20} or as +F<name> defined in an "
              "--mdef file; for frequencies counted (+F) or estimated (+FO) from the data use NONREV";
    else if (params.optimize_from_given_params)
        why = "is not supported with -optfromgiven";
    else if (params.optimize_model_rate_joint)
        why = "is not supported with -jointopt";
    else if (!params.model_joint.empty() || params.link_model)
        why = "does not support --model-joint or --link-model yet";
    else if (params.model_name.find("+T") != string::npos || params.optimize_params_use_hmm)
        why = "is not supported in tree mixtures or HMM models";
    else if (params.model_name.find("+H") != string::npos)
        why = "is not supported with heterotachy (+H) models";
    else if (params.site_freq_file || params.tree_freq_file)
        why = "is not supported with site-specific frequency models";
    else if (params.optimize_linked_gtr)
        why = "is not supported with --link-exchange-rates";
    else if (params.matrix_exp_technique == MET_EIGEN_DECOMPOSITION)
        why = "is not supported with --eigen";
    else if (params.model_name.substr(0, 4) == "TEST" || params.model_name.substr(0, 2) == "MF")
        why = "cannot be tested by ModelFinder; give it with -m";
    else if (params.alisim_active)
        why = "cannot simulate alignments; simulate from its fitted matrix instead";
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

    const char *init = Params::getInstance().model_name_init;
    if (init) {
        NxsModel *nxs = models->findModel(init);
        if (nxs)
            readParametersString(nxs->description, false);
        else
            readParameters(init, false);
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
    int ndim = getNDim();
    double *h = new double[ndim + 1];
    double fx = targetFunk(x);
    for (int dim = 1; dim <= ndim; dim++) {
        double temp = x[dim];
        x[dim] = temp + fixedeq::scaledStep(temp, upper);
        h[dim] = x[dim] - temp;
        dfx[dim] = targetFunk(x);
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

void ModelNonrevFixedEq::checkFit(double optimizer_score) {
    // decisions 009 and 029: re-score the returned point and report; nothing is restored
    fits++;
    optimizer_lnl = optimizer_score;
    decomposeRateMatrix();
    if (build_failed)
        outWarning("NQC: the optimizer returned coordinates that do not give a valid rate matrix; "
                   "the last valid coordinates are kept");
    phylo_tree->clearAllPartialLH();
    rescored_lnl = phylo_tree->computeLikelihood();
    at_lower = at_upper = outside = 0;
    for (int i = 0; i < accepted.size(); i++) {
        at_lower += (accepted[i] == lower);
        at_upper += (accepted[i] == upper);
        outside += (accepted[i] < lower || accepted[i] > upper);
    }
    Map<const Matrix<double, Dynamic, Dynamic, RowMajor> > q(rate_matrix, num_states, num_states);
    residuals = fixedeq::diagnostics(q, target, nu);
    if (verbose_mode >= VB_MED)
        writeInfo(cout);
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
        if (verbose_mode >= VB_DEBUG)
            cout << "NQC: proposal rejected: " << fixedeq::statusText(status) << endl;
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
    startCheckpoint();
    checkpoint->put("chart", string(CHART));
    checkpoint->put("chart_version", CHART_VERSION);
    checkpoint->put("target", numbers17(target.data(), num_states));
    checkpoint->put("references", integers(ref));
    checkpoint->put("domain", numbers17(domain, 2));
    checkpoint->put("derivative", string(DERIVATIVE));
    checkpoint->put("coordinates", numbers17(accepted.data(), (int)accepted.size()));
    checkpoint->put("fits", fits);
    if (fits > 0)
        checkpoint->put("accepted_score", number17(rescored_lnl));
    endCheckpoint();
    ModelProtein::saveCheckpoint();
}

void ModelNonrevFixedEq::restoreCheckpoint() {
    // The base restore decomposes from the coordinates held now; this restore then rebuilds from
    // the saved ones. Absent keys leave the model as it is, as IQ-TREE's own restores do.
    ModelProtein::restoreCheckpoint();
    double domain[2] = {lower, upper};
    string coordinates, chart, saved_target, saved_ref, saved_domain;
    int chart_version = 0;
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
        checkpoint->get("fits", fits);
        string score;
        if (checkpoint->getString("accepted_score", score))
            rescored_lnl = convert_double(score.c_str());
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
    if (phylo_tree)
        phylo_tree->clearAllPartialLH();
}

void ModelNonrevFixedEq::printMatrix(ostream &out, const double *m) {
    for (int i = 0; i < num_states; i++)
        out << numbers17(m + i * num_states, num_states) << endl;
}

void ModelNonrevFixedEq::writeInfo(ostream &out) {
    if (!constructed) {
        ModelProtein::writeInfo(out);
        return;
    }
    string letters;
    for (int r : ref)
        letters += STATES[r];
    out << "NQC target as given (sum " << number17(target_raw_sum) << "): " << target_text << endl;
    out << "NQC target normalized: " << numbers17(target.data(), num_states) << endl;
    out << "NQC reference destinations: " << letters << endl;
    out << "NQC coordinate domain: [" << number17(lower) << ", " << number17(upper) << "]" << endl;
    if (fits == 0)
        return;
    out << "NQC fit " << fits << ": start lnL " << number17(start_lnl) << ", optimizer lnL "
        << number17(optimizer_lnl) << ", re-scored lnL " << number17(rescored_lnl) << endl;
    if (rescored_lnl < start_lnl)
        out << "NQC fit " << fits << ": the re-scored lnL is below the start by "
            << number17(start_lnl - rescored_lnl) << " (reported, not restored)" << endl;
    out << "NQC coordinates at the lower bound " << at_lower << ", at the upper bound " << at_upper
        << ", outside the domain " << outside << endl;
    out << "NQC residuals: |pi*Q| " << number17(residuals.stationarity) << ", row sums "
        << number17(residuals.row_sums) << ", |mean rate - 1| " << number17(residuals.mean_rate_error)
        << ", nu " << number17(residuals.nu_residual) << ", minimum rate " << number17(residuals.min_rate) << endl;
}
