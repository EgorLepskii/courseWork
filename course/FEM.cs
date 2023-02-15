using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace course
{
    public class FEM
    {
        Func<double, double, double> f;
        List<BoundaryCondition> BoundaryConditions;
        Func<double, double, double> lambda;
        Func<double, double, double> gamma;
        List<Elem> elems;
        List<Edge> edges;
        public List<double> xverts;
        public List<double> yverts;
        public List<double> q;
        Matrix mat;

        public FEM(Func<double, double, double> f, List<BoundaryCondition> boundaryConditions, Func<double, double, double> lambda, Func<double, double, double> gamma)
        {
            this.f = f;
            BoundaryConditions = boundaryConditions;
            this.lambda = lambda;
            this.gamma = gamma;
            elems = new();
            edges = new();
            xverts = new();
            yverts = new();
            ReadMesh();
        }

        private void ReadMesh()
        {
            var str = File.ReadAllLines("txt/Verts.txt");
            int vertscount = str.Count();
            foreach (var item in str)
            {
                var cur = item.Split(' ');
                xverts.Add(double.Parse(cur[0]));
                yverts.Add(double.Parse(cur[1]));
            }
            str = File.ReadAllLines("txt/Edges.txt");
            foreach (var item in str)
            {
                var cur = item.Split(' ');
                edges.Add(new Edge(int.Parse(cur[0]), int.Parse(cur[1])));
            }
            str = File.ReadAllLines("txt/elems.txt");
            foreach (var item in str)
            {
                var cur = item.Split(' ');
                elems.Add(new Elem(int.Parse(cur[0]), int.Parse(cur[1]), int.Parse(cur[2]), int.Parse(cur[3]), int.Parse(cur[4]), int.Parse(cur[5]), vertscount));
            }
        }
        private void GenerateProfile()
        {
            mat = new(xverts.Count + edges.Count());
            mat.di = new List<double>(mat.n);
            mat.ia = new List<int>(mat.n + 1);
            mat.ja = new List<int>();
            mat.al = new List<double>();
            mat.b = new List<double>(mat.n);
            List<SortedSet<int>> list = new List<SortedSet<int>>(mat.n);
            for (int i = 0; i < mat.n; i++)
            {
                list.Add(new SortedSet<int>());
            }
            foreach (var item in elems)
            {
                var dic = item.LocalToGlobal;
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        list[dic[i]].Add(dic[j]);
                    }
                }
            }
            mat.ia.Add(0);
            mat.ia.Add(0);
            mat.di.Add(0);
            mat.b.Add(0);
            for (int i = 1; i < mat.n; i++)
            {
                mat.di.Add(0);
                mat.b.Add(0);
                int count = 0;
                foreach (var item in list[i])
                {
                    if (item < i)
                    {
                        mat.ja.Add(item);
                        mat.al.Add(0);
                        count++;
                    }
                }
                mat.ia.Add(mat.ia[i] + count);
            }
        }
        private void AddLocal(Elem elem)
        {
            var d = CalcDet(elem);
            double[] ax = new double[6]//11 12 13 22 32 33
            {
                (yverts[elem.v2] - yverts[elem.v3]) / d * (yverts[elem.v2] - yverts[elem.v3]) / d,
                (yverts[elem.v3] - yverts[elem.v1]) / d * (yverts[elem.v2] - yverts[elem.v3]) / d,
                (yverts[elem.v1] - yverts[elem.v2]) / d * (yverts[elem.v2] - yverts[elem.v3]) / d,
                (yverts[elem.v3] - yverts[elem.v1]) / d * (yverts[elem.v3] - yverts[elem.v1]) / d,
                (yverts[elem.v1] - yverts[elem.v2]) / d * (yverts[elem.v3] - yverts[elem.v1]) / d,
                (yverts[elem.v1] - yverts[elem.v2]) / d * (yverts[elem.v1] - yverts[elem.v2]) / d
            };
            double[] ay = new double[6]
            {
                (xverts[elem.v3] - xverts[elem.v2]) / d * (xverts[elem.v3] - xverts[elem.v2]) / d,
                (xverts[elem.v1] - xverts[elem.v3]) / d * (xverts[elem.v3] - xverts[elem.v2]) / d,
                (xverts[elem.v2] - xverts[elem.v1]) / d * (xverts[elem.v3] - xverts[elem.v2]) / d,
                (xverts[elem.v1] - xverts[elem.v3]) / d * (xverts[elem.v1] - xverts[elem.v3]) / d,
                (xverts[elem.v2] - xverts[elem.v1]) / d * (xverts[elem.v1] - xverts[elem.v3]) / d,
                (xverts[elem.v2] - xverts[elem.v1]) / d * (xverts[elem.v2] - xverts[elem.v1]) / d
            };
            d = Math.Abs(d);
            double gammaloc = gamma((xverts[elem.v1] + xverts[elem.v2] + xverts[elem.v3]) / 3, (yverts[elem.v1] + yverts[elem.v2] + yverts[elem.v3]) / 3);
            double[] lambdas = new double[3]
            {
                lambda(xverts[elem.v1], yverts[elem.v1]),
                lambda(xverts[elem.v2], yverts[elem.v2]),
                lambda(xverts[elem.v3], yverts[elem.v3])

            };
            double[] fs = new double[6]
            {
                f(xverts[elem.v1], yverts[elem.v1]),
                f(xverts[elem.v2], yverts[elem.v2]),
                f(xverts[elem.v3], yverts[elem.v3]),
                f((xverts[elem.v2] + xverts[elem.v3]) / 2, (yverts[elem.v2] + yverts[elem.v3]) / 2),
                f((xverts[elem.v1] + xverts[elem.v3]) / 2, (yverts[elem.v1] + yverts[elem.v3]) / 2),
                f((xverts[elem.v2] + xverts[elem.v1]) / 2, (yverts[elem.v2] + yverts[elem.v1]) / 2)
            };

            var dic = elem.LocalToGlobal;

            for (int i = 0; i < 6; i++)
            {
                for (int m = 0; m < 6; m++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        mat.di[dic[i]] += (ax[m] + ay[m]) * Matrices.GMatr[i][i][m][k] * d * lambdas[k];
                    }
                    mat.b[dic[i]] += fs[m] * Matrices.MMatr[i][m] * d;
                }
                mat.di[dic[i]] += Matrices.MMatr[i][i] * d * gammaloc;
                for (int j = 0; j < i; j++)
                {
                    int max = dic[i] > dic[j] ? dic[i] : dic[j];
                    int min = dic[i] > dic[j] ? dic[j] : dic[i];
                    int index = mat.ja.BinarySearch(mat.ia[max], mat.ia[max + 1] - mat.ia[max], min, default);
                    for (int m = 0; m < 6; m++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            mat.al[index] += lambdas[k] * (ax[m] + ay[m]) * Matrices.GMatr[i][j][m][k] * d;
                        }
                    }
                    mat.al[index] += gammaloc * Matrices.MMatr[i][j] * d;
                }
            }

        }
        private void AddBoundary1()
        {
            foreach (var item in BoundaryConditions.OfType<BoundaryCondition1>())
            {
                var dic = new Dictionary<int, int> { { 0, edges[item.e].v1 }, { 1, item.e + xverts.Count }, { 2, edges[item.e].v2 } };
                var ug = new double[3] { item.ug(xverts[dic[0]], yverts[dic[0]]),
                                         item.ug((xverts[dic[2]] + xverts[dic[0]]) / 2, (yverts[dic[2]] + yverts[dic[0]]) / 2),
                                         item.ug(xverts[dic[2]], yverts[dic[2]]) };
                for (int m = 0; m < 3; m++)
                {
                    mat.di[dic[m]] = 1;
                    mat.b[dic[m]] = ug[m];
                    for (int k = mat.ia[dic[m]]; k < mat.ia[dic[m] + 1]; k++)
                    {
                        mat.b[mat.ja[k]] -= ug[m] * mat.al[k];
                        mat.al[k] = 0;
                    }
                    for (int i = dic[m] + 1; i < mat.n; i++)
                    {
                        int index = mat.ja.BinarySearch(mat.ia[i], mat.ia[i + 1] - mat.ia[i], dic[m], default);
                        if (index >= 0)
                        {
                            mat.b[i] -= mat.al[index] * ug[m];
                            mat.al[index] = 0;
                        }
                    }
                }
            }
        }
        private void AddBoundary2()
        {
            foreach (var item in BoundaryConditions.OfType<BoundaryCondition2>())
            {
                var dic = new Dictionary<int, int> { { 0, edges[item.e].v1 }, { 1, item.e + xverts.Count }, { 2, edges[item.e].v2 } };
                var thetas = new double[3] { item.theta(xverts[dic[0]], yverts[dic[0]]),
                                         item.theta((xverts[dic[2]] + xverts[dic[0]]) / 2, (yverts[dic[2]] + yverts[dic[0]]) / 2),
                                         item.theta(xverts[dic[2]], yverts[dic[2]]) };
                double len = Math.Sqrt((xverts[dic[2]] - xverts[dic[0]]) * (xverts[dic[2]] - xverts[dic[0]]) +
                                       (yverts[dic[2]] - yverts[dic[0]]) * (yverts[dic[2]] - yverts[dic[0]]));
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        mat.b[dic[i]] += thetas[j] * Matrices.BCMatr[i][j] * len;
                    }
                }
            }

        }

        private void AddBoundary3()
        {
            foreach (var item in BoundaryConditions.OfType<BoundaryCondition3>())
            {
                var dic = new Dictionary<int, int> { { 0, edges[item.e].v1 }, { 1, item.e + xverts.Count }, { 2, edges[item.e].v2 } };
                var ubetta = new double[3] { item.ubetta(xverts[dic[0]], yverts[dic[0]]),
                                         item.ubetta((xverts[dic[2]] + xverts[dic[0]]) / 2, (yverts[dic[2]] + yverts[dic[0]]) / 2),
                                         item.ubetta(xverts[dic[2]], yverts[dic[2]]) };
                double len = Math.Sqrt((xverts[dic[2]] - xverts[dic[0]]) * (xverts[dic[2]] - xverts[dic[0]]) +
                                       (yverts[dic[2]] - yverts[dic[0]]) * (yverts[dic[2]] - yverts[dic[0]]));
                for (int i = 0; i < 3; i++)
                {
                    mat.di[dic[i]] += Matrices.BCMatr[i][i] * item.betta * len;
                    for (int j = 0; j < 3; j++)
                    {
                        mat.b[dic[i]] += ubetta[j] * Matrices.BCMatr[i][j] * len * item.betta;
                    }
                    for (int j = 0; j < i; j++)
                    {
                        int max = dic[i] > dic[j] ? dic[i] : dic[j];
                        int min = dic[i] > dic[j] ? dic[j] : dic[i];
                        int index = mat.ja.BinarySearch(mat.ia[max], mat.ia[max + 1] - mat.ia[max], min, default);
                        mat.al[index] += Matrices.BCMatr[i][j] * item.betta * len;
                    }
                }
            }
        }
        public void Solve(double eps, int maxiter)
        {
            GenerateProfile();
            foreach (var item in elems)
            {
                AddLocal(item);
            }
            AddBoundary2();
            AddBoundary3();
            AddBoundary1();
            q = mat.SolveLosLUPrecond(eps, maxiter);
        }
        public double Getsollution(double x, double y)
        {
            double res = 0;
            bool flag = true;
            int i = 0;
            for (i = 0; i < elems.Count && flag; i++)
            {
                double S23 = Math.Abs((xverts[elems[i].v3] - xverts[elems[i].v2]) * (y - yverts[elems[i].v2]) - (x - xverts[elems[i].v2]) * (yverts[elems[i].v3] - yverts[elems[i].v2]));
                double S31 = Math.Abs((xverts[elems[i].v1] - xverts[elems[i].v3]) * (y - yverts[elems[i].v3]) - (x - xverts[elems[i].v3]) * (yverts[elems[i].v1] - yverts[elems[i].v3]));
                double S12 = Math.Abs((xverts[elems[i].v2] - xverts[elems[i].v1]) * (y - yverts[elems[i].v1]) - (x - xverts[elems[i].v1]) * (yverts[elems[i].v2] - yverts[elems[i].v1]));
                if (Math.Abs(Math.Abs(CalcDet(elems[i])) - (S23 + S31 + S12)) <= 1e-7)
                    flag = false;
            }
            i--;
            if (flag)
                throw new Exception();
            var dic = elems[i].LocalToGlobal;
            var L = XtoLTransform(x, y, elems[i]);
            res += L.L1 * (2 * L.L1 - 1) * q[dic[0]];
            res += L.L2 * (2 * L.L2 - 1) * q[dic[1]];
            res += L.L3 * (2 * L.L3 - 1) * q[dic[2]];
            res += 4 * L.L2 * L.L3 * q[dic[3]];
            res += 4 * L.L1 * L.L3 * q[dic[4]];
            res += 4 * L.L1 * L.L2 * q[dic[5]];
            return res;
        }
        public double CalcDet(Elem el)
        {
            return (xverts[el.v2] - xverts[el.v1]) * (yverts[el.v3] - yverts[el.v1]) - (xverts[el.v3] - xverts[el.v1]) * (yverts[el.v2] - yverts[el.v1]);
        }
        public (double L1, double L2, double L3) XtoLTransform(double x, double z, Elem el)
        {
            return ((xverts[el.v2] * yverts[el.v3] - xverts[el.v3] * yverts[el.v2] + x * (yverts[el.v2] - yverts[el.v3]) + z * (xverts[el.v3] - xverts[el.v2])) / CalcDet(el),
                ((xverts[el.v3] * yverts[el.v1] - xverts[el.v1] * yverts[el.v3] + x * (yverts[el.v3] - yverts[el.v1]) + z * (xverts[el.v1] - xverts[el.v3]))) / CalcDet(el),
                ((xverts[el.v1] * yverts[el.v2] - xverts[el.v2] * yverts[el.v1] + x * (yverts[el.v1] - yverts[el.v2]) + z * (xverts[el.v2] - xverts[el.v1]))) / CalcDet(el));
        }
        public (double r, double z) LtoXTransform(double L1, double L2, double L3, Elem el)
        {
            return (L1 * xverts[el.v1] + L2 * xverts[el.v2] + L3 * xverts[el.v3], L1 * yverts[el.v1] + L2 * yverts[el.v2] + L3 * yverts[el.v3]);
        }
        static class Matrices
        {
            public static double[][] MMatr = new double[6][]
            {
                new double[6] {0.01666666666666665, -0.002777777777777775, -0.002777777777777775, -0.01111111111111111, 0, 0, },
                new double[6] {-0.002777777777777775, 0.01666666666666665, -0.002777777777777775, 0, -0.01111111111111111, 0, },
                new double[6]{ -0.002777777777777775, -0.002777777777777775, 0.01666666666666665, 0, 0, -0.01111111111111111, },
                new double[6]{ -0.01111111111111111, 0, 0, 0.08888888888888889, 0.044444444444444446, 0.044444444444444446, },
                new double[6]{ 0, -0.01111111111111111, 0, 0.044444444444444446, 0.08888888888888889, 0.044444444444444446, },
                new double[6]{ 0, 0, -0.01111111111111111, 0.044444444444444446, 0.044444444444444446, 0.08888888888888889, },
            };

            public static double[][][][] GMatr = new double[6][][][]
{
new double[6][][]
{
new double[6][]
{
new double[3]
{
0.30000000000000004, 0.1, 0.1,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.06666666666666665, -0.06666666666666665, -0.033333333333333326,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.06666666666666665, -0.033333333333333326, -0.06666666666666665,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.1, -0.033333333333333326, -0.06666666666666665,
},
new double[3]
{
0.1, -0.06666666666666665, -0.033333333333333326,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0.1, -0.033333333333333326, -0.06666666666666665,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.46666666666666673, 0.1, 0.1,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0.1, -0.06666666666666665, -0.033333333333333326,
},
new double[3]
{
0.46666666666666673, 0.1, 0.1,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
},
new double[6][][]
{
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.06666666666666665, -0.06666666666666665, -0.033333333333333326,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.1, 0.30000000000000004, 0.1,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.033333333333333326, -0.06666666666666665, -0.06666666666666665,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.033333333333333326, 0.1, -0.06666666666666665,
},
new double[3]
{
0.1, 0.46666666666666673, 0.1,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.033333333333333326, 0.1, -0.06666666666666665,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.06666666666666665, 0.1, -0.033333333333333326,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.1, 0.46666666666666673, 0.1,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.06666666666666665, 0.1, -0.033333333333333326,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
},
new double[6][][]
{
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.06666666666666665, -0.033333333333333326, -0.06666666666666665,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.033333333333333326, -0.06666666666666665, -0.06666666666666665,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.1, 0.1, 0.30000000000000004,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.1, 0.1, 0.46666666666666673,
},
new double[3]
{
-0.033333333333333326, -0.06666666666666665, 0.1,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.1, 0.1, 0.46666666666666673,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.06666666666666665, -0.033333333333333326, 0.1,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.033333333333333326, -0.06666666666666665, 0.1,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.06666666666666665, -0.033333333333333326, 0.1,
},
new double[3]
{
0, 0, 0,
},
},
},
new double[6][][]
{
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.1, -0.033333333333333326, -0.06666666666666665,
},
new double[3]
{
0.1, -0.06666666666666665, -0.033333333333333326,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.033333333333333326, 0.1, -0.06666666666666665,
},
new double[3]
{
0.1, 0.46666666666666673, 0.1,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.1, 0.1, 0.46666666666666673,
},
new double[3]
{
-0.033333333333333326, -0.06666666666666665, 0.1,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.26666666666666666, 0.26666666666666666, 0.8,
},
new double[3]
{
0.26666666666666666, 0.5333333333333333, 0.5333333333333333,
},
new double[3]
{
0.26666666666666666, 0.8, 0.26666666666666666,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.26666666666666666, 0.26666666666666666, 0.8,
},
new double[3]
{
0.13333333333333333, 0.26666666666666666, 0.26666666666666666,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.26666666666666666, 0.13333333333333333, 0.26666666666666666,
},
new double[3]
{
0.26666666666666666, 0.26666666666666666, 0.13333333333333333,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.13333333333333333, 0.26666666666666666, 0.26666666666666666,
},
new double[3]
{
0.26666666666666666, 0.8, 0.26666666666666666,
},
new double[3]
{
0.26666666666666666, 0.13333333333333333, 0.26666666666666666,
},
new double[3]
{
0.26666666666666666, 0.26666666666666666, 0.13333333333333333,
},
new double[3]
{
0, 0, 0,
},
},
},
new double[6][][]
{
new double[6][]
{
new double[3]
{
0.1, -0.033333333333333326, -0.06666666666666665,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.46666666666666673, 0.1, 0.1,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.033333333333333326, 0.1, -0.06666666666666665,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.06666666666666665, 0.1, -0.033333333333333326,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.1, 0.1, 0.46666666666666673,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.06666666666666665, -0.033333333333333326, 0.1,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.26666666666666666, 0.26666666666666666, 0.8,
},
new double[3]
{
0.13333333333333333, 0.26666666666666666, 0.26666666666666666,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.26666666666666666, 0.13333333333333333, 0.26666666666666666,
},
new double[3]
{
0.26666666666666666, 0.26666666666666666, 0.13333333333333333,
},
},
new double[6][]
{
new double[3]
{
0.26666666666666666, 0.26666666666666666, 0.8,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.5333333333333333, 0.26666666666666666, 0.5333333333333333,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.8, 0.26666666666666666, 0.26666666666666666,
},
},
new double[6][]
{
new double[3]
{
0.13333333333333333, 0.26666666666666666, 0.26666666666666666,
},
new double[3]
{
0.26666666666666666, 0.13333333333333333, 0.26666666666666666,
},
new double[3]
{
0.26666666666666666, 0.26666666666666666, 0.13333333333333333,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.8, 0.26666666666666666, 0.26666666666666666,
},
new double[3]
{
0, 0, 0,
},
},
},
new double[6][][]
{
new double[6][]
{
new double[3]
{
0.1, -0.06666666666666665, -0.033333333333333326,
},
new double[3]
{
0.46666666666666673, 0.1, 0.1,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.1, 0.46666666666666673, 0.1,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.06666666666666665, 0.1, -0.033333333333333326,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.033333333333333326, -0.06666666666666665, 0.1,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
-0.06666666666666665, -0.033333333333333326, 0.1,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.13333333333333333, 0.26666666666666666, 0.26666666666666666,
},
new double[3]
{
0.26666666666666666, 0.8, 0.26666666666666666,
},
new double[3]
{
0.26666666666666666, 0.13333333333333333, 0.26666666666666666,
},
new double[3]
{
0.26666666666666666, 0.26666666666666666, 0.13333333333333333,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0.13333333333333333, 0.26666666666666666, 0.26666666666666666,
},
new double[3]
{
0.26666666666666666, 0.13333333333333333, 0.26666666666666666,
},
new double[3]
{
0.26666666666666666, 0.26666666666666666, 0.13333333333333333,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.8, 0.26666666666666666, 0.26666666666666666,
},
new double[3]
{
0, 0, 0,
},
},
new double[6][]
{
new double[3]
{
0.26666666666666666, 0.8, 0.26666666666666666,
},
new double[3]
{
0.5333333333333333, 0.5333333333333333, 0.26666666666666666,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0.8, 0.26666666666666666, 0.26666666666666666,
},
new double[3]
{
0, 0, 0,
},
new double[3]
{
0, 0, 0,
},
},
},
};
            public static double[][] BCMatr = new double[3][]
            {
                new double[3]
                {

                        4.0/30,2.0/30,-1.0/30

                },
                new double[3]
                    {
                        2.0/30,16.0/30,2.0/30
                    },
                new double[3]
                    {
                        -1.0/30,2.0/30,4.0/30
                    }
            };
        }
    }
}