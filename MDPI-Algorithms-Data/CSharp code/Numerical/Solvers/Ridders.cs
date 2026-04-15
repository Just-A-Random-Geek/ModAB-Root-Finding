namespace Proektsoft.Numerical
{
    public static partial class Solver
    {
        // Finds the root of "F(x) = 0" within the interval [x1, x2]
        // with the specified precisions - absolute: aTol and relative: rTol,
        // using Ridder's method:
        // C. Ridders, "A new algorithm for computing a single root of a real continuous function"
        // IEEE Transactions on Circuits and Systems, vol. 26, no. 11, pp. 979-980, November 1979
        // https://doi.org/10.1109/TCS.1979.1084580
        // F(x) must be continuous and sign(F(x1)) ≠ sign(F(x2))

        public static double Ridders(Func<double, double> F, double x1, double x2, 
            double aTol = 1e-14, double rTol = 1e-14)
        {
            if (!Initialize(x1, x2, F, out Node p1, out Node p2))
                return double.NaN;

            var x0 = p1.X;
            for (int i = 1; i <= MaxIterations; ++i)
            {
                Node p3 = new(Node.Mid(p1, p2), F);
                double d = p3.Y * p3.Y - p1.Y * p2.Y;
                if (d <= 0)
                    d = p3.X;
                else
                    d = p3.X + (p3.X - p1.X) * Math.Sign(p1.Y) * p3.Y / Math.Sqrt(d);

                var eps = (aTol + rTol * Math.Abs(d)) / 2.0;
                Node p4 = new(d, F);
                if (p4.Y == 0 || Math.Abs(d - x0) <= eps)
                {
                    EvaluationCount = 2 * i + 2;
                    return d;
                }
                x0 = d;
                if (Math.Sign(p1.Y) == Math.Sign(p4.Y))
                {
                    p1 = p4;
                    if (Math.Sign(p2.Y) == Math.Sign(p3.Y))
                        p2 = p3;
                }
                else
                {
                    p2 = p4;
                    if (Math.Sign(p1.Y) == Math.Sign(p3.Y))
                        p1 = p3;
                }
            }
            EvaluationCount = 2 * MaxIterations + 2;
            return double.NaN;
        }

    }
}
