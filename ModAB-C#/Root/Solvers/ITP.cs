namespace Proektsoft.Root
{
    public static partial class Solver
    {
        // Finds the root of "F(x) = 0" within the interval [x1, x2]
        // with the specified precisions - absolute: aTol and relative: rTol,
        // using the ITP method:
        // Oliveira, I. F. D. and Ricardo H. C. Takahashi.
        // “An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality.”
        // ACM Transactions on Mathematical Software (TOMS) 47 (2021): 1 – 24.
        // https://doi.org/10.1145/3423597
        // F(x) must be continuous and sign(F(x1)) ≠ sign(F(x2))

        public static double ITP(Func<double, double> F, double x1, double x2, 
            double aTol = 1e-14, double rTol = 1e-14)
        {
            if (!Initialize(x1, x2, F,out Node p1, out Node p2))
                return double.NaN;

            double span = p2.X - p1.X;
            int n0 = 1;
            double k1 = 0.2 / span, k2 = 2d;
            double eps2 = aTol + rTol * span;
            int nb = (int)Math.Ceiling(Math.Log2(span / eps2));
            int nmax = nb + n0;
            double xb = Node.Mid(p1, p2);
            for (int i = 1; i <= MaxIterations; ++i)
            {
                double xf = Node.Sec(p1, p2);
                double σ = Math.Sign(xb - xf);
                span = p2.X - p1.X;
                double δ = Math.Min(k1 * Math.Pow(span, k2), Math.Abs(xb - xf));
                double xt = xf + σ * δ;
                double rho = Math.Max(0, Math.Min(eps2 * Math.Pow(2d, nmax - i) - span / 2d, Math.Abs(xt - xb)));
                xt = xb - σ * rho;
                Node p = new(xt, F);
                if (Math.Sign(p1.Y) == Math.Sign(p.Y))
                    p1 = p;
                else
                    p2 = p;

                xb = Node.Mid(p1, p2);
                double eps = aTol + rTol * Math.Abs(xb);
                if (p.Y == 0 || p2.X - p1.X < eps)
                {
                    EvaluationCount = i + 2;
                    return xb;
                }
            }
            EvaluationCount = MaxIterations + 2;
            return double.NaN;
        }
    }
}
