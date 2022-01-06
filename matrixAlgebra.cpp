#include <iostream>
#include <cstring>
#include <iomanip>
#include <math.h>
#include <cstdio>
using namespace std;

#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif


template<typename T>
class Matrix
{
protected:
    int rows, columns;
    T **matrix;
public:
    bool cursed = false;
    //Matrix() = delete;

    Matrix(int r, int c){
        this->rows = r;
        this->columns = c;
        this->matrix = new T*[this->rows];				// allocate the rows

        for(int i = 0; i < this->rows; ++i)
        {
            this->matrix[i] = new T[this->columns];		// allocate the columns
        }

        this->nullify();							// set to zero matrix
    }

     ~Matrix(){
        for(int r=0; r<rows; r++){
            delete [] matrix[r];
        }
         delete [] matrix;
         matrix = NULL;
    }

    Matrix(const Matrix<T> &other){
        this->rows = other.rows;
        this->columns = other.columns;
        this->matrix = new T*[this->rows];				// allocate the rows

        for(int i = 0; i < this->rows; ++i)
        {
            this->matrix[i] = new T[this->columns];		// allocate the columns
        }

        for(int i=0; i<this->rows; i++){
            for(int j=0; j<this->columns; j++){
                this->matrix[i][j] = other.matrix[i][j];
            }
        }
    }

    int getNumberOfRows(){
        return this->rows;
    }

    int getNumberOfColumns(){
        return this->columns;
    }
    
    T getOnIndex(int n, int m){
        if(n>=0 && n<rows && m>=0 && m<columns){
            return matrix[n][m];
        }
        return 0;
    }
    
    void setOnIndex(int n, int m, T value){
        matrix[n][m] = value;
    }

    Matrix<T> transposedMatrix(){
        Matrix<T> newMatrix(this->columns, this->rows);
        for(int i=0; i<this->columns; i++){
            for(int j=0; j<this->rows; j++){
                newMatrix.setOnIndex(i, j, this->getOnIndex(j, i));
            }
        }

        return newMatrix;
    }

    void nullify(){
        for(int i=0; i<this->rows; i++){
            for(int j=0; j<this->columns; j++){
                this->matrix[i][j] = 0;
            }
        }
    }

    void input(){
        for(int i=0; i<rows; i++){
            for(int j=0; j<columns; j++){
                cin>>matrix[i][j];
            }
        }
    }

    void output(){
        for(int i=0; i<rows; i++){
            for(int j=0; j<columns; j++){
                cout<<matrix[i][j]<<" ";
            }
            cout<<"\n";
        }
    }

    void outputAugmentedMatrix(const Matrix<T> &other){
        if(this->rows != other.rows) return;

        for(int i=0; i<rows; i++){
            for(int j=0; j<columns; j++){
                cout<<matrix[i][j]<<" ";
            }
            for(int j=0; j<other.columns; j++){
                cout<<other.matrix[i][j]<<" ";
            }
            cout<<"\n";
        }
    }

    Matrix<T>& operator=(Matrix<T> other){
        cursed = other.cursed;
        if(this->cursed){
            return *this;
        }
        if(this->rows != other.getNumberOfRows() || this->columns != other.getNumberOfColumns()){
            cout<<"Error: the dimensional problem occurred\n";
            this->cursed = true;
            return *this;
        }
        for(int i=0; i<rows; i++){
            for(int j=0; j<columns; j++){
                this->setOnIndex(i, j, other.getOnIndex(i, j));
            }
        }
        return *this;
    }

    bool operator==(Matrix<T> other){
        if(this->rows != other.getNumberOfRows() || this->columns != other.getNumberOfColumns()){
            return false;
        }
        for(int i=0; i<this->rows; i++){
            for(int j=0; j<this->columns; j++){
                if(matrix[i][j] != other.getOnIndex(i, j)){
                    return false;
                }
            }
        }

        return true;
    }

    Matrix<T> operator+(Matrix<T>& other){
        Matrix<T> newMatrix(this->rows, this->columns);
        
        if(this->rows != other.getNumberOfRows() || this->columns != other.getNumberOfColumns()){
            cout<<"Error: the dimensional problem occurred\n";
            newMatrix.cursed = true;
            return newMatrix;
        }

        for(int i=0; i<this->rows; i++){
            for(int j=0; j<this->columns; j++){
                T value = this->getOnIndex(i, j) + other.getOnIndex(i, j);
                newMatrix.setOnIndex(i, j, value);
            }
        }

        return newMatrix;
    }

    Matrix<T> operator-(Matrix<T>& other){
        Matrix<T> newMatrix(this->rows, this->columns);
        
        if(this->rows != other.getNumberOfRows() || this->columns != other.getNumberOfColumns()){
            cout<<"Error: the dimensional problem occurred\n";
            newMatrix.cursed = true;
            return newMatrix;
        }

        for(int i=0; i<this->rows; i++){
            for(int j=0; j<this->columns; j++){
                T value = this->getOnIndex(i, j) - other.getOnIndex(i, j);
                newMatrix.setOnIndex(i, j, value);
            }
        }

        return newMatrix;
    }
//TODO: redo this
    Matrix<T> operator*(Matrix<T>& other){
        int otherRows = other.getNumberOfRows();
        int otherColumns = other.getNumberOfColumns();
        Matrix<T> newMatrix(this->rows, otherColumns);
        
        if(columns != otherRows){
            cout<<"Error: the dimensional problem occurred\n";
            newMatrix.cursed = true;
            return newMatrix;
        }
        
        for(int i=0; i<rows; i++){
            for(int j=0; j<otherColumns; j++){
                T sum = T();
                for(int k=0; k<columns; k++){
                    sum += this->getOnIndex(i, k) * other.getOnIndex(k, j);
                }
                newMatrix.setOnIndex(i, j, sum);
            }
        }

        return newMatrix;
    }
};

template<typename T>
class SquareMatrix: public Matrix<T>{
public:
    SquareMatrix(int c): Matrix<T>(c, c){
    }

    T** getPointerToMatrix()
    {
        return this->matrix;
    }

    ~SquareMatrix(){
    }

    void transformToIdentity(){
        for(int i=0; i<this->columns; i++){
            for(int j=0; j<this->rows; j++){
                if(i == j) this->matrix[i][j] = 1;
                else this->matrix[i][j] = T();
            }
        }
    }

    SquareMatrix<T> eliminationMatrix(int n, int m){
        SquareMatrix<T> newMatrix(this->columns);
        if(n > this->columns || m > this->rows){
            cout<<"Error: the dimensional problem occurred\n";
            newMatrix.cursed = true;
            return newMatrix;
        }
        newMatrix.transformToIdentity();
        T e;
        e = -this->getOnIndex(n, m);
        e /= this->getOnIndex(m, m);
        newMatrix.setOnIndex(n, m, e);

        return newMatrix;
    }
    
    SquareMatrix<T> permutationMatrix(int n, int m){
        SquareMatrix<T> newMatrix(this->columns);
        if(n > this->columns || m > this->rows){
            cout<<"Error: the dimensional problem occurred\n";
            newMatrix.cursed = true;
            return newMatrix;
        }
        newMatrix.transformToIdentity();

        newMatrix.setOnIndex(n, n, T());
        newMatrix.setOnIndex(n, m, 1);
        newMatrix.setOnIndex(m, m, T());
        newMatrix.setOnIndex(m, n, 1);

        return newMatrix;
    }
    

    void eliminateDownwards(){
        int k = 1, currentRow = 0;
        while(true){
            /// PERMUTATION
            int absMax=abs(this->matrix[currentRow][currentRow]);
            int idx=currentRow;
            for(int i=currentRow+1; i<this->rows; i++){
                if(absMax < abs(this->matrix[i][currentRow])){
                    absMax = abs(this->matrix[i][currentRow]);
                    idx = i;
                }
            }
            if(idx != currentRow){
                cout<<"step #"<<k++<<": permutation\n";
                *this = this->permutationMatrix(currentRow, idx) * (*this);
                this->output();
            }
            
            /// ELIMINATION
            for(int i=currentRow+1; i<this->rows; i++){
                bool nextIsFine = true;
                for(int j=0; j<=currentRow; j++){
                    if(this->matrix[i][j]){
                        nextIsFine = false;
                        break;
                    }
                }
                if(!nextIsFine){
                    cout<<"step #"<<k++<<": elimination\n";
                    *this = this->eliminationMatrix(i, currentRow) * (*this);
                    this->output();
                }
            }
            
            /// CHECK THE MATRIX
            currentRow++;
            if(currentRow == this->rows){
                break;
            }
        }
    }

    void eliminateDownwardsAugmented(SquareMatrix<T> &other, int &k){
        SquareMatrix<T> newMatrix(this->rows);
        int currentRow = 0;
        while(true){
            /// PERMUTATION
            int absMax=abs(this->matrix[currentRow][currentRow]);
            int idx=currentRow;
            for(int i=currentRow+1; i<this->rows; i++){
                if(absMax < abs(this->matrix[i][currentRow])){
                    absMax = abs(this->matrix[i][currentRow]);
                    idx = i;
                }
            }
            if(idx != currentRow){
                //cout<<"step #"<<k++<<": permutation\n";
                //cout<<"step #"<<k++<<": elimination\n";
                newMatrix = this->permutationMatrix(currentRow, idx);
                other = newMatrix * other;
                *this = newMatrix * (*this);
                //this->outputAugmentedMatrix(other);
            }
            
            /// ELIMINATION
            for(int i=currentRow+1; i<this->rows; i++){
                bool nextIsFine = true;
                for(int j=0; j<=currentRow; j++){
                    if(this->matrix[i][j]){
                        nextIsFine = false;
                        break;
                    }
                }
                if(!nextIsFine){
                    //cout<<"step #"<<k++<<": elimination\n";
                    newMatrix = this->eliminationMatrix(i, currentRow);
                    other = newMatrix * other;
                    *this = newMatrix * (*this);
                    //this->outputAugmentedMatrix(other);
                    //this->output();
                }
            }
            
            /// CHECK THE MATRIX
            currentRow++;
            if(currentRow == this->rows){
                break;
            }
        }
    }

    void eliminateUpwardsAugmented(SquareMatrix<T> &other, int &k){
        SquareMatrix<T> newMatrix(this->rows);
        /// ELIMINATION
        for(int j=this->columns-1; j>0; j--){
            for(int i=j-1; i>=0; i--){
                if(this->matrix[i][j]){
                    //cout<<"step #"<<k++<<": elimination\n";
                    newMatrix = this->eliminationMatrix(i, j);
                    other = newMatrix * other;
                    *this = newMatrix * (*this);
                    //this->outputAugmentedMatrix(other);
                }
            }
        }
    }

    Matrix<T> invertMatrix(){
        SquareMatrix<T> newMatrix(this->rows);
        newMatrix.transformToIdentity();

        //cout<<"step #0: Augmented Matrix\n";
        //this->outputAugmentedMatrix(newMatrix);
        
        int cnt=1;
        //cout<<"Direct way:\n";
        this->eliminateDownwardsAugmented(newMatrix, cnt);
        //cout<<"Way back:\n";
        this->eliminateUpwardsAugmented(newMatrix, cnt);

        //cout<<"Diagonal normalization:\n";
        for(int i=0; i<this->rows; i++){
            for(int j=0; j<this->columns; j++){
                T newValue = newMatrix.getOnIndex(i, j)/this->getOnIndex(i, i);
                newMatrix.setOnIndex(i, j, newValue);
            }
        }
        this->transformToIdentity();
        //this->outputAugmentedMatrix(newMatrix);

        return newMatrix;
    }

    Matrix<T> JacobiMethod(Matrix<T> &B, double &epsilon)
    {
        Matrix<T> res(this->rows, 1);
        for(int i=0; i<this->rows; i++){
            double sum = 0;
            for(int j=0; j<this->rows; j++){
                if(i == j) continue;
                sum += this->getOnIndex(j, i);
            }
            if(sum > this->getOnIndex(i, i)){
                cout<<"The method is not applicable!";
                return res;
            }
        }

        Matrix<T> D(this->rows, 1);
        SquareMatrix<T> A(this->rows);

        for(int i=0; i<this->rows; i++){
            D.setOnIndex(i, 0, 1/this->getOnIndex(i, i));
        }
        
        for(int i=0; i<this->rows; i++){
            for(int j=0; j<this->columns; j++){
                if(i!=j){
                    A.setOnIndex(i, j, -this->getOnIndex(i,j)*D.getOnIndex(i, 0));
                }
            }
        }
        cout<<"alpha:\n";
        A.output();

        for(int i=0; i<this->rows; i++){
            B.setOnIndex(i, 0, B.getOnIndex(i, 0)*D.getOnIndex(i, 0));
            res.setOnIndex(i, 0, B.getOnIndex(i, 0));
        }
        cout<<"beta:\n";
        B.output();

        int cnt=0;
        Matrix<T> oldResult(this->rows, 1);
        Matrix<T> product(this->rows, 1);
        while(true){
            double curEpsinon = 0;

            cout<<"x("<<cnt++<<"):\n";

            oldResult = res;
            res.output();
            product = A*res;
            res = B + product;
            
            for(int i=0; i<this->rows; i++){
                double dif = abs(oldResult.getOnIndex(i, 0)-res.getOnIndex(i, 0));
                curEpsinon +=  dif*dif; 
            }
            curEpsinon = sqrt(curEpsinon);
            cout<<"e: "<<curEpsinon<<"\n";
            

            if(curEpsinon < epsilon){
                break;
            } 
        }
        cout<<"x("<<cnt++<<"):\n";
        res.output();
        return res;
    }

    //TODO: implement Seidel's Method
    Matrix<T> SeidelMethod(Matrix<T> &B, double &epsilon)
    {
        Matrix<T> res(this->rows, 1);
        for(int i=0; i<this->rows; i++){
            double sum = 0;
            for(int j=0; j<this->rows; j++){
                if(i == j) continue;
                sum += this->getOnIndex(j, i);
            }
            if(sum > this->getOnIndex(i, i)){
                cout<<"The method is not applicable!";
                return res;
            }
        }

        Matrix<T> D(this->rows, 1);
        SquareMatrix<T> A(this->rows);

        for(int i=0; i<this->rows; i++){
            D.setOnIndex(i, 0, 1/this->getOnIndex(i, i));
        }
        
        for(int i=0; i<this->rows; i++){
            for(int j=0; j<this->columns; j++){
                if(i!=j){
                    A.setOnIndex(i, j, -this->getOnIndex(i,j)*D.getOnIndex(i, 0));
                }
            }
        }
        cout<<"alpha:\n";
        A.output();

        for(int i=0; i<this->rows; i++){
            B.setOnIndex(i, 0, B.getOnIndex(i, 0)*D.getOnIndex(i, 0));
            res.setOnIndex(i, 0, B.getOnIndex(i, 0));
        }
        cout<<"beta:\n";
        B.output();

        int cnt=0;
        Matrix<T> oldResult(this->rows, 1);
        Matrix<T> product(this->rows, 1);
        while(true){
            double curEpsinon = 0;

            cout<<"x("<<cnt++<<"):\n";

            oldResult = res;
            res.output();
            product = A*res;
            res = B + product;
            
            for(int i=0; i<this->rows; i++){
                double dif = abs(oldResult.getOnIndex(i, 0)-res.getOnIndex(i, 0));
                curEpsinon +=  dif*dif; 
            }
            curEpsinon = sqrt(curEpsinon);
            cout<<"e: "<<curEpsinon<<"\n";
            

            if(curEpsinon < epsilon){
                break;
            } 
        }
        cout<<"x("<<cnt++<<"):\n";
        res.output();
        return res;
    }
    

    T determinantStepByStep(){
        T result = 1;
        int k = 1, currentRow = 0;
        while(true){
            /// PERMUTATION
            int absMax=abs(this->matrix[currentRow][currentRow]);
            int idx=currentRow;
            for(int i=currentRow+1; i<this->rows; i++){
                if(absMax < abs(this->matrix[i][currentRow])){
                    absMax = abs(this->matrix[i][currentRow]);
                    idx = i;
                }
            }
            if(idx != currentRow){
                cout<<"step #"<<k++<<": permutation\n";
                *this = this->permutationMatrix(currentRow, idx) * (*this);
                this->output();
            }
            
            /// ELIMINATION
            for(int i=currentRow+1; i<this->rows; i++){
                bool nextIsFine = true;
                for(int j=0; j<=currentRow; j++){
                    if(this->matrix[i][j]){
                        nextIsFine = false;
                        break;
                    }
                }
                if(!nextIsFine){
                    cout<<"step #"<<k++<<": elimination\n";
                    *this = this->eliminationMatrix(i, currentRow) * (*this);
                    this->output();
                }
            }
            
            /// CHECK THE MATRIX
            currentRow++;
            if(currentRow == this->rows){
                break;
            }
        }

        for(int i=0; i<this->rows; i++){
            result *= this->matrix[i][i];
        }

        return result;
    }

    T determinant(){
        T result = 1;
        
        this->eliminateDownwards();

        for(int i=0; i<this->rows; i++){
            result *= this->matrix[i][i];
        }

        return result;
    }

    SquareMatrix<T>& operator=(SquareMatrix<T> other){
        this->cursed = other.cursed;
        if(this->cursed){
            return *this;
        }
        if(this->rows != other.getNumberOfRows() || this->columns != other.getNumberOfColumns()){
            cout<<"Error: the dimensional problem occurred\n";
            this->cursed = true;
            return *this;
        }
        for(int i=0; i<this->columns; i++){
            for(int j=0; j<this->rows; j++){
                this->matrix[i][j] = other.getOnIndex(i, j);
            }
        }

        return *this;
    }

    SquareMatrix<T>& operator=(Matrix<T> other){
        this->cursed = other.cursed;
        if(this->cursed){
            return *this;
        }
        if(this->rows != other.getNumberOfRows() || this->columns != other.getNumberOfColumns()
                || other.getNumberOfRows() != other.getNumberOfColumns())
        {
            cout<<"Error: the dimensional problem occurred\n";
            this->cursed = true;
            return *this;
        }
        for(int i=0; i<this->rows; i++){
            for(int j=0; j<this->columns; j++){
                this->matrix[i][j] = other.getOnIndex(i, j);
            }
        }

        return *this;
    }

    bool operator==(SquareMatrix<T> other){
        if(this->rows != other.getNumberOfRows() || this->columns != other.getNumberOfColumns()){
            return false;
        }
        for(int i=0; i<this->rows; i++){
            for(int j=0; j<this->columns; j++){
                if(this->matrix[i][j] != other.getOnIndex(i, j)){
                    return false;
                }
            }
        }

        return true;
    }

    SquareMatrix<T> operator+(SquareMatrix<T>& other){
        SquareMatrix<T> newMatrix(this->rows);
        
        if(this->rows != other.getNumberOfRows() || this->columns != other.getNumberOfColumns()){
            cout<<"Error: the dimensional problem occurred\n";
            newMatrix.cursed = true;
            return newMatrix;
        }

        for(int i=0; i<this->rows; i++){
            for(int j=0; j<this->columns; j++){
                T value = this->getOnIndex(i, j) + other.getOnIndex(i, j);
                newMatrix.setOnIndex(i, j, value);
            }
        }

        return newMatrix;
    }

    SquareMatrix<T> operator-(SquareMatrix<T>& other){
        SquareMatrix<T> newMatrix(this->rows);
        
        if(this->rows != other.getNumberOfRows() || this->columns != other.getNumberOfColumns()){
            cout<<"Error: the dimensional problem occurred\n";
            newMatrix.cursed = true;
            return newMatrix;
        }

        for(int i=0; i<this->rows; i++){
            for(int j=0; j<this->columns; j++){
                T value = this->getOnIndex(i, j) - other.getOnIndex(i, j);
                newMatrix.setOnIndex(i, j, value);
            }
        }

        return newMatrix;
    }

    SquareMatrix<T> operator*(const SquareMatrix<T>& other){
        int otherRows = other.rows;
        int otherColumns = other.columns;
        SquareMatrix<T> newMatrix(this->rows);
        
        if(this->columns != otherRows){
            cout<<"Error: the dimensional problem occurred\n";
            newMatrix.cursed = true;
            return newMatrix;
        }
        
        for(int i=0; i<this->rows; i++){
            for(int j=0; j<otherColumns; j++){
                T sum = T();
                for(int k=0; k<this->columns; k++){
                    sum += this->getOnIndex(i, k) * other.matrix[k][j];
                }
                newMatrix.setOnIndex(i, j, sum);
            }
        }

        return newMatrix;
    }
    Matrix<T> operator*(Matrix<T>& other){
        int otherRows = other.getNumberOfRows();
        int otherColumns = other.getNumberOfColumns();
        Matrix<T> newMatrix(this->rows, otherColumns);
        
        if(this->columns != otherRows){
            cout<<"Error: the dimensional problem occurred\n";
            newMatrix.cursed = true;
            return newMatrix;
        }

        for(int i=0; i<this->rows; i++){
            for(int j=0; j<otherColumns; j++){
                T sum = T();
                for(int k=0; k<this->columns; k++){
                    sum += this->getOnIndex(i, k) * other.getOnIndex(k, j);
                }
                newMatrix.setOnIndex(i, j, sum);
            }
        }

        return newMatrix;
    }
};

void solveLeastSquare()
{
    cout<<fixed<<setprecision(4);

    int m;
    cin>>m;
    Matrix<double> t(m, 1), b(m, 1);
    //int t[105], b[105];
    for(int i=0; i<m; i++){
        double a;
        cin>>a;
        t.setOnIndex(i, 0, a);
        cin>>a;
        b.setOnIndex(i, 0, a);
    }
    int n;
    cin>>n;

    Matrix<double> A(m, n+1);
    for(int i=0; i<m; i++){
        A.setOnIndex(i, 0, 1);
        A.setOnIndex(i, 1, t.getOnIndex(i, 0));
        for(int j=2; j<=n; j++){
            int newValue = t.getOnIndex(i, 0) * A.getOnIndex(i, j-1);
            A.setOnIndex(i, j, newValue);
        }
    }
    cout<<"A:\n";
    A.output();
    
    Matrix<double> AT = A.transposedMatrix();
    SquareMatrix<double> ATA(n+1); 
    cout<<"A_T*A:\n";
    
    ATA = AT * A;
    ATA.output();
    
    cout<<"(A_T*A)^-1:\n";
    ATA = ATA.invertMatrix();
    ATA.output();
    
    cout<<"A_T*b:\n";
    Matrix<double> ATb = AT*b;
    ATb.output();

    cout<<"x~:\n";
    Matrix<double> res(n+1, 1);
    res = ATA * ATb;
    res.output();

#ifdef WIN32
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
#else 
    FILE* pipe = popen(GNUPLOT_NAME, "w");
#endif

    if(pipe != NULL){
        fprintf(pipe, "%s\n", "plot '-' title 'line' with lines, '-' title 'dots'");
        
        double prec = 0.01;
        double xMin=1000000, xMax=-1000000;
        for(int i=0; i<m; i++){
            xMin = min(xMin, t.getOnIndex(i, 0));
            xMax = max(xMax, t.getOnIndex(i, 0));
        }
        for(double x=xMin-7; x<xMax+7; x+=prec){
            // cout<<"&";
            double y=0, a=1;
            for(int i=0; i<=n; i++){
                y+=a*res.getOnIndex(i, 0);
                a*=x;
            }
            fprintf(pipe, "%f\t%f\n", x, y);
        }
        fprintf(pipe, "%s\n", "e");
        
        for(int i=0; i<m; i++){
            fprintf(pipe, "%f\t%f\n", t.getOnIndex(i, 0), b.getOnIndex(i, 0));
        }
        fprintf(pipe, "%s\n", "e");

        fflush(pipe);
    }
    else{
        cout<<"could not open pipe";
    }

#ifdef WIN32
    _pclose(pipe);
#else 
    pclose(pipe);
#endif
}


int main()
{
    solveLeastSquare();
}