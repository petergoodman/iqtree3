/*
 * linkedmodelhooks.h
 *
 *  Created on: 23/10/2026
 *      Author: Peter Goodman
 *
 * The calls through which PartitionModel::optimizeLinkedModel lets a linked substitution model
 * check, differentiate and record its fits (docs/agent/DECISIONS.md, decision 005). A model takes
 * part by inheriting this interface; every other model keeps the unchanged path. It is a header of
 * its own because NQC's header is included by modelmixture.cpp and partitionmodel.h pulls in the
 * super tree, so putting it in either would create an include cycle.
 */

#ifndef LINKEDMODELHOOKS_H_
#define LINKEDMODELHOOKS_H_

#include <vector>

class ModelSubst;
class Optimization;

class LinkedModelHooks
{
public:
    virtual ~LinkedModelHooks() {}

    /**
        called before each linked fit; stops with an error if the models cannot share one fit
        @param linked every partition model under the linked name, in the order
               PartitionModel::targetFunk sums them
    */
    virtual void startLinkedFit(const std::vector<ModelSubst*> &linked) = 0;

    /**
        the derivative the linked fit's optimizer uses in place of Optimization::derivativeFunk
        @param objective the PartitionModel whose targetFunk is minimized
        @return the function value at x
    */
    virtual double linkedDerivative(Optimization &objective, double x[], double dfx[]) = 0;

    /**
        called after each linked fit, once the optimizer's point is in place
        @param linked as for startLinkedFit
        @param optimizer_lnl the log-likelihood the optimizer returned
    */
    virtual void endLinkedFit(const std::vector<ModelSubst*> &linked, double optimizer_lnl) = 0;
};

#endif /* LINKEDMODELHOOKS_H_ */
