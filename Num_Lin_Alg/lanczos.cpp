#include <iostream>
#include <Eigen/Dense>

int main() {
    srand(time(0));
    int n = 100;
    int m = 10;
    std::vector<Eigen::VectorXcd> vecs;
    std::vector<double> Bs;
    std::vector<std::complex<double>> As;
    Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(n,n);
    A = (A+A.adjoint())/(2.0);
    Eigen::VectorXcd v = Eigen::VectorXcd::Random(n);
    v = v / v.norm();

    Eigen::VectorXcd w = A*v;
    std::complex<double> alpha = w.dot(v);
    As.push_back(alpha);
    w = w - alpha*v;
    Eigen::VectorXcd v_next;

    for (int i = 1; i<m; i++) {
        double B = w.norm();
        Bs.push_back(B);
        v_next = w / B;
        w = A*v_next - B*v;
        alpha = w.dot(v_next);
        w = w - alpha*v_next;
        As.push_back(alpha);
        v = v_next;
    }
    
    return 0;
}