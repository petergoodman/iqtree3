/*
 * fixedeqchart.cpp
 *
 *  Created on: 19/10/2026
 *      Author: Peter Goodman
 */

#include "fixedeqchart.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <locale>
#include <sstream>

namespace fixedeq {

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

bool validReferences(const std::vector<int> &ref) {
    const int n = static_cast<int>(ref.size());
    if (n < 3)
        return false;
    for (int i = 0; i < n; i++)
        if (ref[i] < 0 || ref[i] >= n || ref[i] == i)
            return false;
    return true;
}

bool validTarget(const VectorXd &pi) {
    for (int i = 0; i < pi.size(); i++)
        if (!std::isfinite(pi(i)) || !(pi(i) > 0.0))
            return false;
    return true;
}

bool square(const MatrixXd &M, int n) { return M.rows() == n && M.cols() == n; }

// PLAN.md constraint 2: every off-diagonal rate finite and strictly positive
Status checkRates(const MatrixXd &Q) {
    for (int i = 0; i < Q.rows(); i++)
        for (int j = 0; j < Q.cols(); j++) {
            if (j == i)
                continue;
            if (!std::isfinite(Q(i, j)))
                return Status::nonFiniteRate;
            if (!(Q(i, j) > 0.0))
                return Status::nonPositiveRate;
        }
    return Status::ok;
}

// checks the off-diagonal rates, then sets each diagonal entry to minus its row's off-diagonal sum
Status completeDiagonal(MatrixXd &Q) {
    const Status s = checkRates(Q);
    if (s != Status::ok)
        return s;
    for (int i = 0; i < Q.rows(); i++) {
        double sum = 0.0;
        for (int j = 0; j < Q.cols(); j++)
            if (j != i)
                sum += Q(i, j);
        if (!std::isfinite(sum))
            return Status::nonFiniteRate;
        Q(i, i) = -sum;
    }
    return Status::ok;
}

bool separator(char c) {
    return c == ',' || c == '/' || c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r';
}

bool digit(char c) { return c >= '0' && c <= '9'; }

// [+-]?(d+(.d*)?|.d+)([eE][+-]?d+)?, ASCII only
bool decimal(const std::string &t) {
    size_t i = 0;
    const size_t n = t.size();
    if (i < n && (t[i] == '+' || t[i] == '-'))
        i++;
    const size_t whole = i;
    while (i < n && digit(t[i]))
        i++;
    bool digits = i > whole;
    if (i < n && t[i] == '.') {
        const size_t frac = ++i;
        while (i < n && digit(t[i]))
            i++;
        digits = digits || i > frac;
    }
    if (!digits)
        return false;
    if (i < n && (t[i] == 'e' || t[i] == 'E')) {
        i++;
        if (i < n && (t[i] == '+' || t[i] == '-'))
            i++;
        const size_t exponent = i;
        while (i < n && digit(t[i]))
            i++;
        if (i == exponent)
            return false;
    }
    return i == n;
}

std::string str(double v) {
    std::ostringstream s;
    s.imbue(std::locale::classic());
    s.precision(17);
    s << v;
    return s.str();
}

}  // namespace

const char *statusText(Status status) {
    switch (status) {
    case Status::ok:
        return "ok";
    case Status::invalidInput:
        return "inputs of the wrong size, or a reference destination out of range or on the diagonal";
    case Status::nonFiniteCoordinate:
        return "a coordinate is not finite";
    case Status::softmaxUnderflow:
        return "a jump probability underflows to zero";
    case Status::nonPositiveNu:
        return "the jump chain's stationary vector is not finite and positive";
    case Status::nonFiniteRate:
        return "a rate is not finite";
    case Status::nonPositiveRate:
        return "an off-diagonal rate is not strictly positive";
    case Status::invalidTarget:
        return "the target is not a finite, strictly positive vector of the right length";
    }
    return "unknown status";
}

std::vector<int> freeDestinations(int i, const std::vector<int> &ref) {
    std::vector<int> free;
    for (int j = 0; j < static_cast<int>(ref.size()); j++)
        if (j != i && j != ref[i])
            free.push_back(j);
    return free;
}

Status jumpMatrix(const VectorXd &z, const std::vector<int> &ref, MatrixXd &K) {
    const int n = static_cast<int>(ref.size());
    if (!validReferences(ref) || z.size() != n * (n - 2))
        return Status::invalidInput;
    if (!z.allFinite())
        return Status::nonFiniteCoordinate;
    MatrixXd out = MatrixXd::Zero(n, n);
    VectorXd logit(n);
    int k = 0;
    for (int i = 0; i < n; i++) {
        // the reference logit is 0; subtracting the row maximum keeps every exp at most 1
        double top = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i)
                continue;
            logit(j) = j == ref[i] ? 0.0 : z(k++);
            top = std::max(top, logit(j));
        }
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i)
                continue;
            out(i, j) = std::exp(logit(j) - top);
            if (out(i, j) == 0.0)
                return Status::softmaxUnderflow;
            sum += out(i, j);
        }
        // a subnormal exp can still round to zero once divided by a row sum above 1
        for (int j = 0; j < n; j++) {
            if (j == i)
                continue;
            out(i, j) /= sum;
            if (out(i, j) == 0.0)
                return Status::softmaxUnderflow;
        }
    }
    K = out;
    return Status::ok;
}

Status solveStationary(const MatrixXd &M, VectorXd &x) {
    const int n = static_cast<int>(M.rows());
    if (n < 1 || M.cols() != n)
        return Status::invalidInput;
    MatrixXd A(n + 1, n);
    A.row(0).setOnes();
    A.bottomRows(n) = M.transpose();
    VectorXd b = VectorXd::Zero(n + 1);
    b(0) = 1.0;
    const VectorXd solution = A.colPivHouseholderQr().solve(b);
    if (!solution.allFinite())
        return Status::nonPositiveNu;
    x = solution;
    return Status::ok;
}

Status generatorFromJump(const MatrixXd &K, const VectorXd &pi, MatrixXd &Q, VectorXd &nu) {
    const int n = static_cast<int>(pi.size());
    if (n < 3 || !square(K, n) || K.diagonal().cwiseAbs().maxCoeff() != 0.0)
        return Status::invalidInput;
    if (!validTarget(pi))
        return Status::invalidTarget;
    VectorXd v;
    Status s = solveStationary(K - MatrixXd::Identity(n, n), v);
    if (s != Status::ok)
        return s;
    for (int i = 0; i < n; i++)
        if (!(v(i) > 0.0))
            return Status::nonPositiveNu;
    MatrixXd out(n, n);
    for (int i = 0; i < n; i++) {
        const double exit = v(i) / pi(i);
        for (int j = 0; j < n; j++)
            out(i, j) = j == i ? 0.0 : exit * K(i, j);
    }
    s = completeDiagonal(out);
    if (s != Status::ok)
        return s;
    Q = out;
    nu = v;
    return Status::ok;
}

Status build(const VectorXd &z, const VectorXd &pi, const std::vector<int> &ref, MatrixXd &Q, VectorXd &nu) {
    MatrixXd K;
    const Status s = jumpMatrix(z, ref, K);
    if (s != Status::ok)
        return s;
    return generatorFromJump(K, pi, Q, nu);
}

Status inverse(const MatrixXd &Q, const std::vector<int> &ref, VectorXd &z) {
    const int n = static_cast<int>(ref.size());
    if (!validReferences(ref) || !square(Q, n))
        return Status::invalidInput;
    Status s = checkRates(Q);
    if (s != Status::ok)
        return s;
    VectorXd out(n * (n - 2));
    int k = 0;
    for (int i = 0; i < n; i++) {
        const double exit = -Q(i, i);
        if (!std::isfinite(exit))
            return Status::nonFiniteRate;
        if (!(exit > 0.0))
            return Status::nonPositiveRate;
        const double reference = Q(i, ref[i]) / exit;
        for (int j = 0; j < n; j++) {
            if (j == i || j == ref[i])
                continue;
            out(k) = std::log((Q(i, j) / exit) / reference);
            if (!std::isfinite(out(k)))
                return Status::nonFiniteCoordinate;
            k++;
        }
    }
    z = out;
    return Status::ok;
}

std::vector<int> references(const MatrixXd &Q) {
    const int n = static_cast<int>(Q.rows());
    std::vector<int> ref(n, -1);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (j != i && (ref[i] < 0 || Q(i, j) > Q(i, ref[i])))
                ref[i] = j;
    return ref;
}

Status reversibleSeed(const MatrixXd &R, const VectorXd &pi, MatrixXd &Q) {
    const int n = static_cast<int>(pi.size());
    if (n < 3 || !square(R, n))
        return Status::invalidInput;
    if (!validTarget(pi))
        return Status::invalidTarget;
    MatrixXd out(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            out(i, j) = j == i ? 0.0 : R(i, j) * pi(j);
    Status s = completeDiagonal(out);
    if (s != Status::ok)
        return s;
    double rate = 0.0;
    for (int i = 0; i < n; i++)
        rate += pi(i) * -out(i, i);
    out /= rate;
    s = checkRates(out);
    if (s != Status::ok)
        return s;
    Q = out;
    return Status::ok;
}

Status transferSeed(const MatrixXd &Q0, const VectorXd &pi, MatrixXd &Q) {
    const int n = static_cast<int>(pi.size());
    if (n < 3 || !square(Q0, n))
        return Status::invalidInput;
    Status s = checkRates(Q0);
    if (s != Status::ok)
        return s;
    MatrixXd K(n, n);
    for (int i = 0; i < n; i++) {
        const double exit = -Q0(i, i);
        if (!std::isfinite(exit))
            return Status::nonFiniteRate;
        if (!(exit > 0.0))
            return Status::nonPositiveRate;
        for (int j = 0; j < n; j++)
            K(i, j) = j == i ? 0.0 : Q0(i, j) / exit;
    }
    VectorXd nu;
    return generatorFromJump(K, pi, Q, nu);
}

Diagnostics diagnostics(const MatrixXd &Q, const VectorXd &pi, const VectorXd &nu) {
    const int n = static_cast<int>(pi.size());
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Diagnostics d = {nan, nan, nan, nan, nan, nan, nan};
    if (n < 1 || !square(Q, n) || nu.size() != n)
        return d;
    const VectorXd exit = -Q.diagonal();
    d.stationarity = (pi.transpose() * Q).cwiseAbs().maxCoeff();
    d.row_sums = (Q.rowwise().sum().cwiseAbs().array() / exit.array()).maxCoeff();
    d.mean_rate_error = std::fabs(pi.dot(exit) - 1.0);
    d.min_rate = d.min_flux = std::numeric_limits<double>::infinity();
    MatrixXd K = Q;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            if (j != i) {
                d.min_rate = std::min(d.min_rate, Q(i, j));
                d.min_flux = std::min(d.min_flux, pi(i) * Q(i, j));
            }
        K.row(i) /= exit(i);
        K(i, i) = 0.0;
    }
    d.nu_residual = (nu.transpose() * K - nu.transpose()).cwiseAbs().maxCoeff();
    const VectorXd expected = nu.cwiseQuotient(pi);
    d.exit_rate_mismatch = ((exit - expected).cwiseAbs().array() / expected.array()).maxCoeff();
    return d;
}

Status parseTarget(const std::string &text, int count, double min_freq, VectorXd &pi, double &raw_sum,
                   std::string &error) {
    std::vector<std::string> tokens;
    std::string token;
    for (char c : text + ",") {
        if (separator(c)) {
            if (!token.empty())
                tokens.push_back(token);
            token.clear();
        } else {
            token += c;
        }
    }
    if (static_cast<int>(tokens.size()) != count) {
        error = "the target has " + std::to_string(tokens.size()) + " entries; it needs " + std::to_string(count);
        return Status::invalidTarget;
    }
    std::vector<double> values;
    for (size_t k = 0; k < tokens.size(); k++) {
        const std::string where = "entry " + std::to_string(k + 1) + " (" + tokens[k] + ")";
        if (!decimal(tokens[k])) {
            error = where + " is not a decimal number";
            return Status::invalidTarget;
        }
        std::istringstream s(tokens[k]);
        s.imbue(std::locale::classic());
        double v = 0.0;
        s >> v;
        if (s.fail() || s.peek() != std::char_traits<char>::eof() || !std::isfinite(v)) {
            error = where + " is not a finite number";
            return Status::invalidTarget;
        }
        if (!(v > 0.0) || v < min_freq) {
            error = where + " must be greater than 0 and at least " + str(min_freq);
            return Status::invalidTarget;
        }
        values.push_back(v);
    }
    double sum = 0.0;
    for (double v : values)
        sum += v;
    if (!(std::fabs(sum - 1.0) <= 1e-6)) {
        error = "the entries sum to " + str(sum) + ", not within 1e-6 of 1";
        return Status::invalidTarget;
    }
    VectorXd out(count);
    for (int i = 0; i < count; i++)
        out(i) = values[i] / sum;
    pi = out;
    raw_sum = sum;
    return Status::ok;
}

double scaledStep(double x, double upper, double eta, double scale) {
    double h = eta * std::max(scale, std::fabs(x));
    if (x + h > upper)
        h = -h;
    const double probe = x + h;
    return probe - x;
}

}  // namespace fixedeq
