using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace course
{
    public class Matrix
    {
        public int n;
        public List<int> ia = new();
        public List<int> ja = new();
        public List<double> di = new();
        public List<double> al = new();
        public List<double> b = new();
        public List<double> di_LU = new();
        public List<double> au_LU = new();
        public List<double> al_LU = new();
        public Matrix(int n)
        {
            this.n = n;
        }
        public void Clear()
        {
            for (int i = 0; i < n; i++)
            {
                di[i] = 0;
                di_LU[i] = 0;
                b[i] = 0;
            }
            int m = al.Count;
            for (int i = 0; i < m; i++)
            {
                al[i] = 0;
                al_LU[i] = 0;
                au_LU[i] = 0;
            }
        }
        public void LU()
        {
            foreach (var item in di)
            {
                di_LU.Add(item);
            }
            foreach (var item in al)
            {
                au_LU.Add(item);
                al_LU.Add(item);
            }
            for (int i = 0; i < n; i++)
            {
                di_LU[i] = di[i];
            }
            var q = al.Count;
            for (int i = 0; i < q; i++)
            {
                al_LU[i] = al[i];
                au_LU[i] = al[i];
            }
            for (int i = 0; i < n; i++)
            {
                double sumdi = 0.0;

                int i0 = ia[i];
                int i1 = ia[i + 1];


                for (int k = i0; k < i1; k++)
                {
                    int j = ja[k];
                    int j0 = ia[j];
                    int j1 = ia[j + 1];
                    int ik = i0;
                    int kj = j0;

                    double suml = 0.0;
                    double sumu = 0.0;

                    while (ik < k)
                    {

                        if (ja[ik] == ja[kj])
                        {

                            suml += al_LU[ik] * au_LU[kj];
                            sumu += au_LU[ik] * al_LU[kj];
                            ik++;
                            kj++;
                        }

                        else
                        {
                            if (ja[ik] > ja[kj])
                            {
                                kj++;
                            }
                            else
                            {
                                ik++;
                            }
                        }
                    }

                    al_LU[k] = al_LU[k] - suml;
                    au_LU[k] = (au_LU[k] - sumu) / di_LU[j];
                    sumdi += al_LU[k] * au_LU[k];
                }

                di_LU[i] = di_LU[i] - sumdi;
            }
        }
        List<double> LUDirect(List<double> rpart)
        {
            List<double> res = new();
            for (int i = 0; i < n; i++)
            {
                res.Add(rpart[i]);
            }

            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int j = ia[i]; j < ia[i + 1]; j++)
                    sum += al_LU[j] * res[ja[j]];
                res[i] -= sum;
                res[i] /= di_LU[i];
            }
            return res;
        }
        List<double> LUReverse(List<double> rpart)
        {
            List<double> res = new();
            for (int i = 0; i < n; i++)
            {
                res.Add(rpart[i]);
            }

            for (int i = n - 1; i >= 0; i--)
            {
                for (int j = ia[i]; j < ia[i + 1]; j++)
                    res[ja[j]] -= au_LU[j] * res[i];
            }
            return res;
        }
        public List<double> LoS_precond(List<double> x0, double eps, int maxiter)
        {
            int k = 1;
            List<double> buf = MatrixMult(x0);
            double bnorm = 0;
            for (int i = 0; i < n; i++)
            {
                buf[i] = b[i] - buf[i];
            }
            double rnorm = Math.Sqrt(DotProduct(buf, buf));
            List<double> r = LUDirect(buf);
            bnorm = Math.Sqrt(DotProduct(b, b));
            List<double> z = LUReverse(r);
            buf = MatrixMult(z);
            List<double> p = LUDirect(buf);
            double resid = 1;
            while (resid > eps && k < maxiter)
            {
                double pp = DotProduct(p, p);
                double pr = DotProduct(p, r);
                double alpha = pr / pp;
                for (int i = 0; i < n; i++)
                {
                    x0[i] += alpha * z[i];
                    r[i] -= alpha * p[i];
                }
                rnorm = Math.Sqrt(DotProduct(r, r));
                List<double> Ur = LUReverse(r);
                buf = MatrixMult(Ur);
                buf = LUDirect(buf);
                double betta = -(DotProduct(p, buf) / pp);
                for (int i = 0; i < n; i++)
                {
                    z[i] = Ur[i] + betta * z[i];
                    p[i] = buf[i] + betta * p[i];
                }
                double test1 = 0;
                double test2 = 0;
                var asd = MatrixMult(x0);
                for (int i = 0; i < n; i++)
                {
                    test1 += (asd[i] - b[i]) * (asd[i] - b[i]);
                    test2 += b[i] * b[i];
                }
                resid = Math.Sqrt(test1 / test2);
                k++;
            }
            //Console.WriteLine($"{k} {rnorm / bnorm} {resid}");
            return x0;
        }
        public List<double> LOS(List<double> x0,  double eps, int maxiter)
        {
            double bnorm = Math.Sqrt(DotProduct(b, b));

            List<double> Ar = new();
            List<double> r = MatrixMult(x0);
            List<double> z = new();
            for (int i = 0; i < n; i++)
            {
                r[i] = b[i] - r[i];
                z.Add(r[i]);
            }
            List<double> p = MatrixMult(z);
            int k = 0;
            double alpha, betta, rnorm = Math.Sqrt(DotProduct(r, r));
            while (k < maxiter && rnorm / bnorm > eps)
            {
                alpha = DotProduct(p, r) / DotProduct(p, p);
                for (int i = 0; i < n; i++)
                {
                    x0[i] += alpha * z[i];
                    r[i] -= alpha * p[i];
                }
                Ar = MatrixMult(r);
                betta = -DotProduct(p, Ar) / DotProduct(p, p);
                rnorm = Math.Sqrt(DotProduct(r, r));
                for (int i = 0; i < n; i++)
                {
                    z[i] = r[i] + betta * z[i];
                    p[i] = Ar[i] + betta * p[i];
                }
                k++;
            }
            Console.WriteLine($"{rnorm / bnorm} {k}");
            return x0;
        }
        double DotProduct(List<double> vec1, List<double> vec2)
        {
            if (vec1.Count != vec2.Count)
                throw new Exception();
            double res = 0;
            int m = vec1.Count;
            for (int i = 0; i < m; i++)
            {
                res += vec1[i] * vec2[i];
            }
            return res;
        }
        List<double> MatrixMult(List<double> x)
        {
            if (x.Count != n)
                throw new Exception();
            List<double> res = new List<double>(x.Count);
            for (int i = 0; i < n; i++)
            {
                res.Add(0);
            }
            for (int i = 0; i < n; i++)
            {
                res[i] = x[i] * di[i];
                for (int k = ia[i]; k < ia[i + 1]; k++)
                {
                    int j = ja[k];
                    res[i] += al[k] * x[j];
                    res[j] += al[k] * x[i];
                }
            }
            return res;
        }
        public List<double> SolveLosLUPrecond(double eps,int maxiter)
        {
            List<double> x0 = new();
            for (int i = 0; i < n; i++)
            {
                x0.Add(0);
            }
            LU();
            return LoS_precond(x0, eps, maxiter);
        }
    }
}

