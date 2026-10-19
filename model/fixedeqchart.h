/*
 * fixedeqchart.h
 *
 *  Created on: 19/10/2026
 *      Author: Peter Goodman
 *
 * The jump-chain log-ratio chart of the fixed-equilibrium non-reversible model (D01 in
 * docs/agent/DECISIONS.md). For a strictly positive target pi*, it maps n(n-2) unconstrained
 * coordinates one-to-one onto the generators Q with positive off-diagonal rates, zero row sums,
 * pi* Q = 0 and unit mean rate at pi*:
 *
 *   z -> K   row softmax over the destinations j != i, with the reference logit z_{i,r(i)} = 0
 *   K -> nu  nu K = nu, sum(nu) = 1
 *   q_ij = (nu_i / pi*_i) K_ij,  q_ii = -sum_{j != i} q_ij
 *
 * Coordinates are row-major; within row i they run over the destinations other than i and r(i)
 * in increasing index (freeDestinations). Pure functions of their arguments, with no IQ-TREE
 * headers, no static mutable state and no process exit, so that they can be tested alone
 * (unittest/) and called inside PartitionModel's OpenMP loop. A function that can fail returns a
 * Status and writes its outputs only when it returns Status::ok.
 */

#ifndef FIXEDEQCHART_H_
#define FIXEDEQCHART_H_

#include <string>
#include <vector>

#include <Eigen/Dense>

namespace fixedeq {

enum class Status {
    ok,
    invalidInput,          // wrong sizes, or a reference that is out of range or on the diagonal
    nonFiniteCoordinate,
    softmaxUnderflow,      // a jump probability rounds to exactly zero
    nonPositiveNu,         // the jump chain's stationary vector is not finite and positive
    nonFiniteRate,
    nonPositiveRate,       // an off-diagonal rate is not strictly positive
    invalidTarget
};

const char *statusText(Status status);

/** derivative step of D02 and decision 005: eta max(scale, |x|) */
const double STEP_ETA = 1e-4;
const double STEP_SCALE = 1.0;

/** the destinations of row i other than i and ref[i], increasing: the order of row i's coordinates */
std::vector<int> freeDestinations(int i, const std::vector<int> &ref);

/** the jump matrix K of coordinates z (size n(n-2)), with zero diagonal */
Status jumpMatrix(const Eigen::VectorXd &z, const std::vector<int> &ref, Eigen::MatrixXd &K);

/** x with x M = 0 and sum(x) = 1, from [1^T; M^T] x = e_1 by column-pivoted Householder QR, as
    computeStateFreqFromQMatrix solves pi from Q (model/modelmarkov.cpp) */
Status solveStationary(const Eigen::MatrixXd &M, Eigen::VectorXd &x);

/** Q and nu of a jump matrix K at target pi */
Status generatorFromJump(const Eigen::MatrixXd &K, const Eigen::VectorXd &pi, Eigen::MatrixXd &Q,
                         Eigen::VectorXd &nu);

/** jumpMatrix followed by generatorFromJump */
Status build(const Eigen::VectorXd &z, const Eigen::VectorXd &pi, const std::vector<int> &ref,
             Eigen::MatrixXd &Q, Eigen::VectorXd &nu);

/** the coordinates of Q: K_ij = q_ij / (-q_ii), z_ij = log(K_ij / K_{i,r(i)}); no box is applied */
Status inverse(const Eigen::MatrixXd &Q, const std::vector<int> &ref, Eigen::VectorXd &z);

/** the row maxima of Q's off-diagonal entries, ties to the lower index (decision 011) */
std::vector<int> references(const Eigen::MatrixXd &Q);

/** q_ij = R_ij pi_j for symmetric exchangeabilities R, normalized to unit mean rate (D05) */
Status reversibleSeed(const Eigen::MatrixXd &R, const Eigen::VectorXd &pi, Eigen::MatrixXd &Q);

/** Q0's jump chain rebuilt at pi: the flux-preserving transfer of a non-reversible seed */
Status transferSeed(const Eigen::MatrixXd &Q0, const Eigen::VectorXd &pi, Eigen::MatrixXd &Q);

/** residuals of a generator against target pi and the nu it was built with; reported, never used
    to reject (decision 026) */
struct Diagnostics {
    double stationarity;        // max_j |(pi Q)_j|
    double row_sums;            // max_i |sum_j q_ij| / (-q_ii)
    double mean_rate_error;     // |sum_i pi_i (-q_ii) - 1|
    double min_rate;            // min off-diagonal q_ij
    double min_flux;            // min off-diagonal pi_i q_ij
    double nu_residual;         // max_j |(nu K)_j - nu_j|, K = q_ij / (-q_ii)
    double exit_rate_mismatch;  // max_i |-q_ii - nu_i / pi_i| / (nu_i / pi_i)
};

Diagnostics diagnostics(const Eigen::MatrixXd &Q, const Eigen::VectorXd &pi, const Eigen::VectorXd &nu);

/**
    Read a target written as `count` decimal numbers separated by commas, ASCII whitespace or
    slashes (decision 006). Every entry must be greater than 0 and at least min_freq, and the sum,
    taken left to right, within 1e-6 of 1; the entries are then divided by that sum once. On
    success writes pi and the raw sum; otherwise writes only error.
*/
Status parseTarget(const std::string &text, int count, double min_freq, Eigen::VectorXd &pi,
                   double &raw_sum, std::string &error);

/**
    Signed forward-difference step for coordinate x (D02, decision 005): h0 = eta max(scale, |x|),
    negated when x + h0 > upper, returned as the displacement (x + h0) - x that the probe actually
    makes, as Optimization::derivativeFunk does. That displacement is exactly h0 when |x| >= h0.
*/
double scaledStep(double x, double upper, double eta = STEP_ETA, double scale = STEP_SCALE);

}  // namespace fixedeq

#endif /* FIXEDEQCHART_H_ */
