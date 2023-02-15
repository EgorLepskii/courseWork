namespace course
{
    internal class Program
    {
        static void Main(string[] args)
        {
            List<BoundaryCondition> BCs = new();

            BCs.Add(new BoundaryCondition1(1, (x, y) => 1));
            BCs.Add(new BoundaryCondition3(3, 1, (x, y) => 21));
            BCs.Add(new BoundaryCondition2(0, (x, y) => 0));
            BCs.Add(new BoundaryCondition2(4, (x, y) => 0));

            FEM fem = new FEM((x, y) => -4 + 4 * x * x, BCs, (x, y) => 2, (x, y) => 4);
            fem.Solve(3e-15, 10000);

            Console.WriteLine($"{fem.Getsollution(1, 1)} {fem.Getsollution(2, 1)} {fem.Getsollution(3, 1)}");

            /*for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    Console.WriteLine($"{1 + 1.0 / 8 + i / 2.0} {1.0 / 8 + j / 2.0} {fem.Getsollution(1 + 1.0 / 8 + i / 2.0, 1.0 / 8 + j / 2.0)}");
                    Console.WriteLine($"{1 + 3.0 / 8 + i / 2.0} {3.0 / 8 + j / 2.0} {fem.Getsollution(1 + 3.0 / 8 + i / 2.0, 3.0 / 8 + j / 2.0)}");
                }
            }*/
            Console.WriteLine("Hello World!");
        }
    }
}