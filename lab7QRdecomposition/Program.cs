using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
class Program
{
	static void Main(string[] args)
	{
		Matrix<double> A = DenseMatrix.OfArray(new double[,]
		{
				{ 0.42016, -16937, 0.10087, -0.28570 },
				{ 0.19439, -0.76571, 0.45605, -0.13218 },
				{ -0.61729, 0.28952, -0.17253, 0.41974 },
				{ -0.20038, 0.78932, -0.47011, 0.13625 }
		});

		Console.WriteLine("Исходная матрица:\n");
		for(int i = 0; i < A.RowCount; i++)
        {
			for(int j = 0; j < A.ColumnCount; j++)
            {
				Console.Write("\t\t" + A.ToArray()[i, j] + "  ");
            }
			
			Console.WriteLine();
        }

		var qr = A.QR();

		Console.WriteLine("Q:\n");
		for (int i = 0; i < qr.Q.RowCount; i++)
		{
			for (int j = 0; j < qr.Q.ColumnCount; j++)
			{
				Console.Write("\t" + qr.Q.ToArray()[i, j] + "  ");
			}

			Console.WriteLine();
		}


		Console.WriteLine("R:\n");
		for (int i = 0; i < qr.R.RowCount; i++)
		{
			for (int j = 0; j < qr.R.ColumnCount; j++)
			{
				Console.Write("\t\t{0:0.#}  ", qr.R.ToArray()[i, j]);
			}

			Console.WriteLine();
		}

		double[,] y = new double[4, 5];

		y[0, 0] = 1;
		Console.WriteLine("\nВекторы\n");
		for (int k = 0; k < 4; k++)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					y[i, k + 1] += A[i, j] * y[j, k];
				}
			}
		}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 5; j++)
			{
				Console.Write("\t\t{0:0.0000}", y[i, j]);
			}

			Console.WriteLine();
		}


		for (int i = 0; i < 4; i++)
		{
			double tmp;
			tmp = y[i, i];
			for (int j = 4; j >= i; j--)
            {
				y[i, j] /= tmp;
			}

			for (int j = i + 1; j < 4; j++)
			{
				tmp = y[j, i];
				for (int k = 4; k >= i; k--)
                {
					y[j, k] -= tmp * y[i, k];
				}
			}
		}

		double[] xx = new double[4];
		double[] p = new double[5];

		xx[3] = y[3, 4];
		for (int i = 2; i >= 0; i--)
		{
			xx[i] = y[i, 4];
			for (int j = i + 1; j < 4; j++)
			{
				xx[i] -= y[i, j] * xx[j];
			}

		}

		for (int i = 4; i > 0; i--)
		{
			p[i] = xx[4 - i];
		}

		p[0] = 1;

		double r;
		double[] a = new double[4];
		double[] b = new double[4]; 
		
		for (int i = 0; i < 4; i++)
		{ 
			r = 0;
			for (int j = 0; j < 4; j++)
			{
				r += Math.Abs(A[i, j]);
			}

			r -= A[i, i];
			a[i] = A[i, i] - r;
			b[i] = A[i, i] + r;
		}


		double[] x = new double[4];
		double eps = 0.0001;
		Console.WriteLine("\nСобственные значения матрицы:\n");

		for (int i = 0; i < 4; i++)
		{
			x[i] = Secants(a[i], b[i], eps, p);
			Console.Write("\t\t{0:0.0000}", x[i]);
		}

		Console.ReadKey();
	}


	static double Function(double[] p, double y)
	{
		return p[0] * Math.Pow(y, 4) - p[1] * Math.Pow(y, 3) - p[2] * Math.Pow(y, 2) - p[3] * y - p[4];
	}

	static double Secants(double a, double b, double eps, double[] p)
	{
		while (Math.Abs(b - a) > eps)
		{
			a = b - (b - a) * Function(p, b) / (Function(p, b) - Function(p, a));
			b = a - (a - b) * Function(p, a) / (Function(p, a) - Function(p, b));
		}

		return (a + b) / 2;
	}
}