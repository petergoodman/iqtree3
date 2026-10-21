/*
 * modelnonrevfixedeq.h
 *
 *  Created on: 20/10/2026
 *      Author: Peter Goodman
 *
 * NQC: the non-reversible amino-acid model whose stationary distribution is fixed to a target
 * pi* given by the user (docs/agent/PLAN.md). Its 360 free parameters are the jump-chain
 * log-ratio coordinates of fixedeqchart.h, and every Q it decomposes is built from them, so
 * pi* Q = 0 at every likelihood evaluation and pi* is the root distribution. pi* is never
 * re-solved from Q or changed after construction.
 */

#ifndef MODELNONREVFIXEDEQ_H_
#define MODELNONREVFIXEDEQ_H_

#include <vector>

#include <Eigen/Dense>

#include "fixedeqchart.h"
#include "modelprotein.h"

class ModelNonrevFixedEq : public ModelProtein
{
public:
    /** the public model name, provisional (decision 002) */
    static const char *const NAME;

    /** @return true if model_name is NQC's name, in any case */
    static bool isName(const string &model_name);

    /** stop with an error if model is an NQC model; called where mixture components are created
        (decision 032) */
    static void refuseAsMixtureComponent(ModelSubst *model);

    /**
        constructor (decisions 002, 006, 011, 031)
        @param model_params text in braces after the name; must be empty
        @param freq frequency type from the model string; must be user-defined
        @param freq_params the target pi*, as written in +F{...} or defined for +F<name>
        @param tree associated phylogenetic tree
        @param models_block the built-in and user-defined models
    */
    ModelNonrevFixedEq(string model_params, StateFreqType freq, string freq_params, PhyloTree *tree,
                       ModelsBlock *models_block);

    /** -lnL at x, or 1e30 when the coordinates do not build a valid Q (PLAN.md constraint 2) */
    virtual double targetFunk(double x[]);

    /** forward differences with the scaled step of decision 005 */
    virtual double derivativeFunk(double x[], double dfx[]);

    virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

    /** ModelMarkov::optimizeParameters followed by the post-fit checks of decisions 009 and 029 */
    virtual double optimizeParameters(double gradient_epsilon);

    /** build Q from the coordinates, then the unchanged base decomposition (decision 015) */
    virtual void decomposeRateMatrixNonrev();

    /** pi*, without the base class's renormalization (decision 034) */
    virtual void getStateFrequency(double *state_freq, int mixture = 0);

    /** these accept only pi* itself */
    virtual void setStateFrequency(double *state_freq);
    virtual void adaptStateFrequency(double *state_freq);

    /** these stop with an error, since every caller is refused (decision 034) */
    virtual void setRateMatrix(double *rate_mat);
    virtual void setFullRateMatrix(double *rate_mat, double *freq);

    virtual void startCheckpoint();
    virtual void saveCheckpoint();
    virtual void restoreCheckpoint();

    virtual void writeInfo(ostream &out);

protected:
    virtual void setVariables(double *variables);
    virtual bool getVariables(double *variables);

private:
    /** false until the constructor finishes; until then the overrides defer to the base class,
        because IQ-TREE's model readers decompose and print while NQC reads its seed */
    bool constructed;

    ModelsBlock *models;
    string target_text;
    Eigen::VectorXd target;
    double target_raw_sum;
    std::vector<int> ref;
    double lower, upper;

    /** coordinates last received from the optimizer, and the last ones that built a valid Q,
        which are the published state (CODE_PLAN.md section 2.2) */
    Eigen::VectorXd trial, accepted;
    bool build_failed;
    Eigen::MatrixXd Q;
    Eigen::VectorXd nu;

    /** the post-fit record */
    int fits;
    bool fit_starting;
    bool gradient_printed;
    double start_lnl, optimizer_lnl, rescored_lnl;
    int at_lower, at_upper, outside;
    fixedeq::Diagnostics residuals;

    void refuseUnsupported(const string &model_params, StateFreqType freq, const string &freq_params);
    void readSeed(Eigen::MatrixXd &seed);
    void checkFit(double optimizer_score);
    void printMatrix(ostream &out, const double *m);
};

#endif /* MODELNONREVFIXEDEQ_H_ */
