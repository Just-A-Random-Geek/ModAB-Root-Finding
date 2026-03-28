namespace Proektsoft.Root
{
    public static partial class Solver
    {
        // Finds the root of "F(x) = 0" within the interval [x1, x2]
        // with the specified precisions - absolute: aTol and relative: rTol,
        // using the bisection method
        // F(x) must be continuous and sign(F(x1)) ≠ sign(F(x2))

        public static double Bisection(Func<double, double> F, double x1, double x2,
            double aTol = 1e-14, double rTol = 1e-14)
        {
            if (!Initialize(x1, x2, F, out Node p1, out Node p2))
                return double.NaN;

            for (int i = 1; i <= MaxIterations; ++i)
            {
                Node p3 = new(Node.Mid(p1, p2), F);
                double d = p2.X - p1.X;
                var eps2 = aTol + rTol * Math.Abs(p3.X);
                if (p3.Y == 0.0 || d <= eps2)
                {
                    if (d <= eps2)
                        p3.X = Node.Sec(p1, p2);

                    EvaluationCount = i + 2;
                    return p3.X;
                }
                if (Math.Sign(p1.Y) == Math.Sign(p3.Y))
                    p1 = p3;
                else
                    p2 = p3;
            }
            EvaluationCount = MaxIterations + 2;
            return double.NaN;
        }
    }
}
