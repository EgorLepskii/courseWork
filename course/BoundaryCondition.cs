using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace course
{
    public class BoundaryCondition
    {
        public BoundaryCondition(int edge)
        {
            this.e = edge;
        }

        public int e { get; init; }
    }

    public class BoundaryCondition1 : BoundaryCondition
    {
        public Func<double, double, double> ug;

        public BoundaryCondition1(int edge, Func<double, double, double> ug) : base(edge)
        {
            this.ug = ug;
        }
    }

    public class BoundaryCondition2 : BoundaryCondition
    {
        public Func<double, double, double> theta;

        public BoundaryCondition2(int edge, Func<double, double, double> theta) : base(edge)
        {
            this.theta = theta;
        }
    }
    public class BoundaryCondition3 : BoundaryCondition
    {
        public double betta;
        public Func<double, double, double> ubetta;

        public BoundaryCondition3(int edge, double betta, Func<double, double, double> ubetta) : base(edge)
        {
            this.betta = betta;
            this.ubetta = ubetta;
        }
    }

}