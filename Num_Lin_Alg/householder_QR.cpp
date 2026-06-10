#include <iostream>
#include <Eigen/Dense>

double sign(double x) {
    return (x >= 0) ? 1.0 : -1.0; 
}

Eigen::VectorXd Householder(Eigen::VectorXd a) {
    double norm = a.norm();
    Eigen::VectorXd v = a;
    if (norm != 0) {
        v(0) = a(0) + sign(a(0))*norm;
    }
    return v;
}

Eigen::MatrixXd HPA(Eigen::MatrixXd A, Eigen::VectorXd v) {
    double Beta = 2.0 / v.dot(v);
    Eigen::VectorXd w = Beta*(A.transpose()*v);
    A = A - v*(w.transpose());
    return A;
}

std::pair<Eigen::MatrixXd, std::vector<double>> HQR(Eigen::MatrixXd A) {
    int n = A.cols();
    Eigen::VectorXd v = Eigen::VectorXd::Zero(n);
    std::vector<double> r;
    int m = A.col(0).size();
    for (int i = 0; i<n; i++) {
        v.segment(i,m-i) = Householder(A.col(i).segment(i,m-i));
        A.bottomRightCorner(m-i,n-i) = HPA(A.bottomRightCorner(m-i,n-i),v.segment(i,m-i));
        r.push_back(A(i,i));
        A.col(i).segment(i,m-i) = v.segment(i,m-i);
    }
    return {A,r};
}

Eigen::MatrixXd find_R(Eigen::MatrixXd A, std::vector<double> r) {
    Eigen::VectorXd p = Eigen::Map<Eigen::VectorXd>(r.data(),r.size());
    Eigen::MatrixXd Out = p.asDiagonal();
    for (int i = 0; i < A.col(0).size(); i++) {
        Out.col(i).segment(0,i) = A.col(i).segment(0,i);
    }
    return Out;
}

Eigen::MatrixXd find_Q(Eigen::MatrixXd A) {
    int n = A.row(0).size();
    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(n,n);
    Eigen::VectorXd v = A.col(0);
    Eigen::MatrixXd Temp = v*v.transpose();
    Q = Q - (2.0/(v.dot(v)))*Temp;
    for (int i = 1; i < n; i++) {
        v = Eigen::VectorXd::Zero(n);
        v.segment(i,n-i) = A.col(i).segment(i,n-i);
        Temp = Eigen::MatrixXd::Identity(n,n) - (2.0/(v.dot(v)))*(v*v.transpose());
        Q = Q*Temp;
    }
    return Q;

}

int main() {
    int n = 5;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
    Eigen::MatrixXd A_original = A;
    std::vector<double> r;
    std::tie(A, r) = HQR(A);
    Eigen::Matrix Q = find_Q(A);
    Eigen::MatrixXd R = find_R(A,r);
    std::cout << "Q:" << std::endl;
    std::cout << Q << std::endl;
    std::cout << "R:" << std::endl;
    std::cout << R << std::endl;
    std::cout << "QQ^T" << std::endl;
    std::cout << Q*Q.transpose()<< std::endl;
    std::cout << "QR - A" << std::endl;
    std::cout << Q*R - A_original;
}
