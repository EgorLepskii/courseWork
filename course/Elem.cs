using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace course
{
    public class Elem
    {
        public Elem(int v1, int v2, int v3, int e1, int e2, int e3, int countOfVerts)
        {
            this.v1 = v1;
            this.v2 = v2;
            this.v3 = v3;
            this.e1 = e1;
            this.e2 = e2;
            this.e3 = e3;
            LocalToGlobal = new Dictionary<int, int>() { { 0, v1 }, { 1, v2 }, { 2, v3 }, { 3, e1 + countOfVerts }, { 4, e2 + countOfVerts }, { 5, e3 + countOfVerts } };
        }

        public IReadOnlyDictionary<int, int> LocalToGlobal { get; init; }
        public int v1 { get; init; }
        public int v2 { get; init; }
        public int v3 { get; init; }
        public int e1 { get; init; }
        public int e2 { get; init; }
        public int e3 { get; init; }
    }
}