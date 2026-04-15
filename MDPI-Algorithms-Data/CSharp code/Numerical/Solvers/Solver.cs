namespace Proektsoft.Numerical
{
    // This class provides methods for solving the equation F(x) = 0 numerically
    // All methods are of bracketing type. They require an initial interval to be specified
    // that contains the root. The function F(x) must be continuous within the interval
    // and sign(F(x1)) ≠ sign(F(x2))

    public static partial class Solver
    {

        private const int MaxIterations = 200;

        public static int EvaluationCount { get; private set; }


        private static bool Initialize(double x1, double x2, Func<double, double> F,
            out Node p1, out Node p2)
        {
            if (x1 > x2)
            {
                p1 = new(x2, F);
                p2 = new(x1, F);
            }
            else
            {
                p1 = new(x1, F);
                p2 = new(x2, F);
            }
            EvaluationCount = 0;
            return Math.Sign(p1.Y) != Math.Sign(p2.Y);
        }
    }
}