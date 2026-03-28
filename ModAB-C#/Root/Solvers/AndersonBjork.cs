namespace Proektsoft.Root
{
    public static partial class Solver
    {
        // Finds the root of "F(x) = 0" within the interval [x1, x2]
        // with the specified precisions - absolute: aTol and relative: rTol,
        // using Anderson-Bjork's method:
        // Anderson, N., Björck, Å. A new high order method of regula-falsi type
        // for computing a root of an equation.
        // BIT 13, 253–264 (1973). 
        // https://doi.org/10.1007/BF01951936
        // F(x) must be continuous and sign(F(x1)) ≠ sign(F(x2))

        public static double AndersonBjork(Func<double, double> F, double x1, double x2,
            double aTol = 1e-14, double rTol = 1e-14)
        {
            if (!Initialize(x1, x2, F, out Node p1, out Node p2))
                return double.NaN;

            int side = 0;
            var x0 = p1.X;
            for (int i = 1; i <= MaxIterations; ++i)
            {
                Node p3 = new(Node.Sec(p1, p2), F);
                var eps = (aTol + rTol * Math.Abs(p3.X)) / 2.0;
                if (p3.Y == 0 || Math.Abs(p3.X - x0) <= eps)
                {
                    EvaluationCount = i + 2;
                    return p3.X;
                }
                x0 = p3.X;
                if (Math.Sign(p1.Y) == Math.Sign(p3.Y))
                {
                    if (side == 1)
                    {
                        double m = 1 - p3.Y / p1.Y;
                        if (m <= 0)
                            p2.Y /= 2;
                        else
                            p2.Y *= m;
                    }
                    else
                        side = 1;
                    
                    p1 = p3;
                }
                else
                {
                    if (side == -1)
                    {
                        double m = 1 - p3.Y / p2.Y;
                        if (m <= 0)
                            p1.Y /= 2;
                        else
                            p1.Y *= m;
                    }
                    else
                        side = -1;
                    
                    p2 = p3;
                }
            }
            EvaluationCount = MaxIterations + 2;
            return double.NaN;
        }
    }
}
