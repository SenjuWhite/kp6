
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra.Factorization;
using static System.Net.Mime.MediaTypeNames;
double[,] test1 = new double[,] {
    {1,2,3},
    {4,2,1},
    {5,2,2},
};
double[,] A15 = new double[,] {
    { 0.5,1.2,2,1},
    { 1.2,2,0.5,1.2},
    { 2,0.5,1,0.5},
    { 1,1.2,0.5,1.6}
};
double[,] A16 = new double[,]{
    {1.8,1.6,1.7,1.8},
    {1.6,2.8,1.5,1.3},
    {1.7,1.5,3.8,1.4},
    {1.8,1.3,1.4,4.8}
};
double[,] A23 = new double[,]
{
{1.6,1.6,1.7,1.8},
    {1.6,2.6,1.3,1.3},
    {1.7,1.3,3.6,1.4},
    {1.8,1.3,1.4,4.6}
};
double[,] A24 = new double[,]
{
    {2,1.6,1.7,1.8},
    {1.6,3,1.7,1.3},
    {1.7,1.7,4,1.4},
    {1.8,1.3,1.4,5}
};
Matrix<double> A_24 = Matrix<double>.Build.DenseOfArray(A24);
Matrix<double> A_23 = Matrix<double>.Build.DenseOfArray(A23);
Matrix<double> A_16 = Matrix<double>.Build.DenseOfArray(A16);
Matrix<double> A_15 = Matrix<double>.Build.DenseOfArray(A15);
void Task2_1(double[,] matrix)
{

    Matrix<double> L = Matrix<double>.Build.DenseOfArray(LU(matrix,4,"l"));
    Matrix<double> U = Matrix<double>.Build.DenseOfArray(LU(matrix, 4, "u"));
    int count = 0;
    var A = L.Multiply(U);
    int iterations = 1000;
    double norm_prev = 0;
    for (int i = 0; i < iterations; i++)
    {      
        var A_new = Matrix<double>.Build.DenseOfArray(LU(A.ToArray(), 4, "u")).Multiply(Matrix<double>.Build.DenseOfArray(LU(A.ToArray(), 4, "l")));      
        double norm = A_new.L2Norm();
        if (Math.Abs(norm - norm_prev) < 0.001)
        {
            break;
        }
        A = A_new;
        norm_prev = norm;
        count++;
    }
    Console.WriteLine($"Кількість ітерацій:{count}");
    Console.WriteLine(A.ToString());
}
void Task2_2(Matrix<double> matrix)
{
    int count = 0;
    var A = matrix;
    int iterations = 1000;
    double norm_prev = double.MaxValue;
    for (int i = 0; i < iterations; i++)
    {
        var qr = A.QR();
        var A_new = qr.R.Multiply(qr.Q);
        double norm = A_new.QR().R.L2Norm();
        if (Math.Abs(A[0, 0] - A_new[0,0])< 0.001&& Math.Abs(A[1, 1] - A_new[1, 1]) < 0.001&& Math.Abs(A[2, 2] - A_new[2, 2]) < 0.001&& Math.Abs(A[3, 3] - A_new[3, 3]) < 0.001)
        {
            break;
        }
        A = A_new;
        norm_prev = norm;
        count++;
    }
    Console.WriteLine($"Кількість ітерацій:{count}");
    Console.WriteLine(A.ToString());
    Console.WriteLine(A.QR().R.ToString());
    Console.WriteLine(A.QR().Q.ToString());
}
Task2_2(A_16);
Task2_1(A16);
double[,] LU(double[,] matrix, int n, string L_or_U )
    {
        double[,] lower = new double[n, n];
        double[,] upper = new double[n, n];

        // Decomposing matrix into Upper and Lower
        // triangular matrix
        for (int i = 0; i < n; i++)
        {
            // Upper Triangular
            for (int k = i; k < n; k++)
            {
                // Summation of L(i, j) * U(j, k)
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (lower[i, j] * upper[j, k]);

                // Evaluating U(i, k)
                upper[i, k] = matrix[i, k] - sum;
            }

            // Lower Triangular
            for (int k = i; k < n; k++)
            {
                if (i == k)
                    lower[i, i] = 1; // Diagonal as 1
                else
                {
                    // Summation of L(k, j) * U(j, i)
                    double sum = 0;
                    for (int j = 0; j < i; j++)
                        sum += (lower[k, j] * upper[j, i]);

                    // Evaluating L(k, i)
                    lower[k, i]
                        = (matrix[k, i] - sum) / upper[i, i];
                }
            }
        }
        return L_or_U == "l" ? lower : upper;
    }


