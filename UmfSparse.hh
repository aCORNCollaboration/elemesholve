#include <map>
using std::map;
#include <vector>
using std::vector;
using std::pair;
#include <stdio.h>

/// Convenience interface for umfpack sparse matrices
class UmfSparse {
public:
    /// Constructor for m x n matrix (dynamically resizeable)
    UmfSparse(int m = 0, int n = 0): mmax(m), nmax(n) { }
    
    /// Destructor
    ~UmfSparse();
    
    /// Element access operator
    double& operator()(int i, int j);
    /// Const access operator
    double operator()(int i, int j) const;
    
    /// Sort data by row
    void sort();
    /// Prepare solver
    void setupSolver();
    /// Solve Ax = b. Resizes x as needed.
    void solve(vector<double>& x, const vector<double>& b);
    
    /// Multiply b = Ax; resizes x as needed.
    void mul(vector<double>& b, const vector<double>& x);
    
    void display() const { printf("Sparse (%i x %i) matrix with %zu entries.\n", mmax, nmax, Ai.size()); }
    
protected:
    vector<int> Ai;             ///< row indices
    vector<int> Aj;             ///< column indices
    vector<double> Ax;          ///< value at row/column
    
    vector<int> Ap;             ///< sorted row data pointers
    bool is_sorted = false;     ///< whether data arrays are in sorted order
    
    map< pair<int, int>, int > index;
    int mmax, nmax;
    void* Symbolic = NULL;
    void* Numeric = NULL;
};

