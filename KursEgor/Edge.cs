using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace KursEgor
{
    public class Edge
    {
        public Edge(int v1, int v2)
        {
            this.v1 = v1;
            this.v2 = v2;
        }

        public int v1 { get; init; }
        public int v2 { get; init; }
    }
}
