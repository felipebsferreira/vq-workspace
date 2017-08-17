// See "Enhanced Fireworks Algorithm", S. Zheng, A. Janecek, Y. Tan (2013)
using System;
using System.Collections.Generic;

namespace FireworksAlgorithm
{
  class FireworksProgram
  {
    static void Main(string[] args)
    {
      Console.WriteLine("\nBegin fireworks algorithm optimization demo\n");
      Console.WriteLine("Goal is to find solution to Ackley's function for 10 variables");
      Console.WriteLine("Function has known min value = 0.0 at (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)");

      int dim = 10; // often number weights to find in ML training
      int n = 5;    // number fireworks
      int maxEpochs = 1000;
      Console.WriteLine("\nSetting Ackley's dimension to " + dim);
      Console.WriteLine("Setting number fireworks to " + n);
      Console.WriteLine("Setting maxEpochs to " + maxEpochs);

      Console.WriteLine("\nBegin algorithm\n");
      double[] bestPosition = Solve(dim, n, maxEpochs); // use fireworks algorithm
      Console.WriteLine("\nAlgorithm complete");

      Console.WriteLine("\nBest solution found: ");
      for (int i = 0; i < dim; ++i)
        Console.Write(bestPosition[i].ToString("F3") + " ");
      Console.WriteLine();

      double error = Error(bestPosition);
      Console.WriteLine("\nError of best solution found = " + error.ToString("F5"));

      Console.WriteLine("\nEnd fireworks algorithm optimization demo\n");
      Console.ReadLine();
    } // Main

    public static double Error(double[] position)
    {
      // Ackley's function has min = 0.0 at (0.0, 0.0, . . 0.0)
      double ss = 0.0;
      for (int i = 0; i < position.Length; ++i)
        ss = ss + (position[i] * position[i]);
      double mss = ss / position.Length;

      double sc = 0.0;
      for (int i = 0; i < position.Length; ++i)
        sc = sc + Math.Cos(2.0 * Math.PI * position[i]);
      double msc = sc / position.Length;

      double f = -20.0 * Math.Exp(-0.2 * Math.Sqrt(mss)) - Math.Exp(msc) + 20.0 + Math.E;
      double trueMinVal = 0.0;
      //return Math.Abs(f - trueMinVal);
      return (f - trueMinVal) * (f - trueMinVal);
    }

    public static double[] Solve(int dim, int n, int maxEpochs)
    {
      int m = n * 10;  // total number of regular sparks for all fireworks
      int mHat = 5;    // number of Gaussian sparks
      double a = 0.04; // controls min sparks per firework
      double b = 0.8;  // controls max sparks per firework
      int A = 40;      // max amplitude
      double minX = -10.0;
      double maxX = 10.0;

      // generate n random fireworks
      Random rnd = new Random(3); // seed = 3 gives representative demo
      Info[] fireworks = new Info[n];
      for (int i = 0; i < n; ++i)
      {
        fireworks[i] = new Info();
        fireworks[i].position = new double[dim];
        for (int j = 0; j < dim; ++j)
          fireworks[i].position[j] = (maxX - minX) * rnd.NextDouble() + minX;
        fireworks[i].error = Error(fireworks[i].position);
      }
 
      // initialize best position, error to dummy values
      double[] bestPosition = new double[dim];
      for (int k = 0; k < dim; ++k)
        bestPosition[k] = fireworks[0].position[k];
      double bestError = fireworks[0].error; // arbitrary

      List<Info>[] sparksList = new List<Info>[n]; // each cell is a ref to a List of Info objects.
      for (int i = 0; i < n; ++i)
        sparksList[i] = new List<Info>(); // empty List

      // main processing loop
      int epoch = 0;
      while (epoch < maxEpochs)
      {
        if (epoch % 100 == 0) // show progress every 100 iterations
        {
          Console.Write("epoch = " + epoch);
          Console.WriteLine(" error at best known position = " + bestError.ToString("F4"));
        }

        int[] numberSparks = NumberSparks(fireworks, m, a, b); // number regular sparks for each firework
        double[] amplitudes = Amplitudes(fireworks, A, epoch, maxEpochs, minX, maxX); // amplitude each firework. epoch, maxEpochs, minX, maxX to establish a min amplitude

        for (int i = 0; i < n; ++i)
          sparksList[i].Clear(); // number of sparks changes each epoch

        for (int i = 0; i < n; ++i) // generate regular sparks for each firework
        {
          double amp = amplitudes[i]; // amplitude for curr firework
          int ns = numberSparks[i];   // number sparks for curr firework

          for (int j = 0; j < ns; ++j) // each spark for curr firework
          {
            Info spark = new Info(); // a spark has a position and error
            spark.position = new double[dim]; // allocate space (ctor doesn't)
            for (int k = 0; k < dim; ++k) // spark position based on its parent firework
              spark.position[k] = fireworks[i].position[k];
            int z = (int)Math.Round(dim * rnd.NextDouble()); // number of random dimensions
            int[] dimensions = PickDimensions(dim, z, epoch); // select z random dimensions. epoch is a rnd seed
            for (int ii = 0; ii < dimensions.Length; ++ii) // each of the randomly selected dimensions
            {
              double h = amp * 2 * rnd.NextDouble() - 1; // displacement hi = +1, lo = -1, (hi - lo) * r + lo
              int k = dimensions[ii]; // convenience
              spark.position[k] += h; // displace from parent firework
              if (spark.position[k] < minX || spark.position[k] > maxX) // bring out-of-range values back in
                spark.position[k] = (maxX - minX) * rnd.NextDouble() + minX;
            }
            spark.error = Error(spark.position);
            sparksList[i].Add(spark);

            // is curr spark global best?
            if (spark.error < bestError)
            {
              bestError = spark.error;
              for (int k = 0; k < dim; ++k)
                bestPosition[k] = spark.position[k];
            }
          } // each new regular spark
        } // each firework

        // pretty hideous parameter passing . . 
        AddGaussianSparks(fireworks, sparksList, dim, mHat, epoch, minX, maxX, bestPosition, ref bestError, rnd);

        // now pick numFireworks new fireworks from all sparks
        // use best spark, worst spark, and n-2 random sparks
        // find best and worst spark
        double[] bestSparkPos = new double[dim];
        double bestSparkErr = double.MaxValue; 

        double[] worstSparkPos = new double[dim];
        double worstSparkErr = double.MinValue;

        for (int i = 0; i < n; ++i) // numFireworks
        {
          for (int j = 0; j < sparksList[i].Count; ++j) // number sparks in each firework
          {
            if (sparksList[i][j].error < bestSparkErr)
            {
              bestSparkErr = sparksList[i][j].error;
              for (int k = 0; k < sparksList[i][j].position.Length; ++k)
                bestSparkPos[k] = sparksList[i][j].position[k];
            }
            if (sparksList[i][j].error > worstSparkErr)
            {
              worstSparkErr = sparksList[i][j].error;
              for (int k = 0; k < sparksList[i][j].position.Length; ++k)
                worstSparkPos[k] = sparksList[i][j].position[k];
            }
          } // each spark
        } // each firework

        for (int k = 0; k < dim; ++k) // first new firework is best spark
          fireworks[0].position[k] = bestSparkPos[k];
        fireworks[0].error = bestSparkErr;

        for (int k = 0; k < dim; ++k) // second new firework is worst spark
           fireworks[1].position[k] = worstSparkPos[k];
        fireworks[1].error = worstSparkErr;

        for (int i= 2; i < n; ++i) // n-2 random sparks 
        {
          int row = rnd.Next(0, n); // more likely to be a good spark
          int cols = sparksList[row].Count;
          int col = rnd.Next(0, cols);
          for (int k = 0; k < dim; ++k)
            fireworks[i].position[k] = sparksList[row][col].position[k];
          fireworks[i].error = sparksList[row][col].error;
        }
        // consider, instead selecting a row based on row-count . . 

        ++epoch;
      } // main loop
      return bestPosition;
    } // Solve

    //private static void ShowVector(int[] vector)
    //{
    //  Console.WriteLine("");
    //  for (int i = 0; i < vector.Length; ++i)
    //    Console.Write(vector[i] + " ");
    //  Console.WriteLine("\n");
    //}

    //private static void ShowInfo(Info info)
    //{
    //  Console.WriteLine("iiiiiiiiiiiiiiiiiiiiiiiiiiiiii");
    //  Console.Write(" pos = ");
    //  for (int k = 0; k < info.position.Length; ++k)
    //    Console.Write(info.position[k].ToString("F2") + " ");
    //  Console.WriteLine("  error = " + info.error.ToString("F3"));

    //  Console.WriteLine("iiiiiiiiiiiiiiiiiiiiiiiiiiiiii");
    //}

    //private static void ShowFireworks(Info[] fireworks)
    //{
    //  Console.WriteLine("ffffffffffffffffffffffffffffff");
    //  for (int i = 0; i < fireworks.Length; ++i)
    //  {
    //    Console.Write(i + " pos = ");
    //    for (int k = 0; k < fireworks[i].position.Length; ++k)
    //      Console.Write(fireworks[i].position[k].ToString("F2") + " ");
    //    Console.WriteLine("  error = " + fireworks[i].error.ToString("F3"));
    //  }
    //  Console.WriteLine("ffffffffffffffffffffffffffffff");
    //}

    //private static void ShowSparksList(List<Info>[] sparks)
    //{
    //  Console.WriteLine("++++++++++++++++++++++++++++");
    //  for (int i = 0; i < sparks.Length; ++i)
    //  {
    //    Console.Write(i + " ");
    //    for (int j = 0; j < sparks[i].Count; ++j)
    //    {
    //      Console.Write(" { ");
    //      for (int k = 0; k < sparks[i][j].position.Length; ++k)
    //        Console.Write(sparks[i][j].position[k].ToString("F2") + " ");
    //      Console.Write(" } ");
    //      Console.WriteLine(" err = " + sparks[i][j].error.ToString("F3"));
    //    }
    //  }
    //  Console.WriteLine("++++++++++++++++++++++++++++");
    //}

    private static int[] PickDimensions(int dim, int z, int seed)
    {
      // pick z random dimensions of a position
      int[] result = new int[z];
      int[] indices = new int[dim];
      for (int i = 0; i < dim; ++i)
        indices[i] = i;

      Random rnd = new Random(seed); // shuffle indices
      for (int i = 0; i < indices.Length; ++i)
      {
        int ri = rnd.Next(i, indices.Length);
        int tmp = indices[ri];
        indices[ri] = indices[i];
        indices[i] = tmp;
      }
      // copy first z indices to result
      for (int i = 0; i < z; ++i)
        result[i] = indices[i];

      return result;
    }

    private static double YMax(Info[] fireworks)
    {
      // largest (worst) error in any firework
      double result = fireworks[0].error;
      for (int i = 1; i < fireworks.Length; ++i)
        if (fireworks[i].error > result)
          result = fireworks[i].error;
      return result;
    }

    private static double YMin(Info[] fireworks)
    {
      // smallest (best) error in any firework
      double result = fireworks[0].error;
      for (int i = 1; i < fireworks.Length; ++i)
        if (fireworks[i].error < result)
          result = fireworks[i].error;
      return result;
    }

    private static int[] NumberSparks(Info[] fireworks, int m, double a, double b)
    {
      // number sparks for each firework
      int n = fireworks.Length;
      int minSparks = (int)Math.Round(a * m); // if n=5, m=50, a=.04, -> 2
      if (minSparks < 1) minSparks = 1;
      int maxSparks = (int)Math.Round(b * m); // if n=5, m=50, b=.8 -> 40
      if (maxSparks > m - (n - 1) * minSparks)
        maxSparks = m - (n - 1) * minSparks;

      double yMax = YMax(fireworks);
      double sumDeltas = 0.0; // sum diffs between yMax and each error
      for (int i = 0; i < n; ++i)
        sumDeltas += yMax - fireworks[i].error;

      int[] numSparks = new int[n]; // the result
      for (int i = 0; i < n; ++i)
      {
        numSparks[i] = (int)Math.Round(m * (yMax - fireworks[i].error + 1.0E-10) / (sumDeltas + 1.0E-10));
        if (numSparks[i] < minSparks)
          numSparks[i] = minSparks;
        else if (numSparks[i] > maxSparks)
          numSparks[i] = maxSparks;
      }
      return numSparks;
    }

    private static double[] Amplitudes(Info[] fireworks, int A, int epoch, int maxEpochs, double minX, double maxX)
    {
      int n = fireworks.Length;
      double yMin = YMin(fireworks);
      double sumDeltas = 0.0; // sum  diffs between yMin and each error
      for (int i = 0; i < n; ++i)
        sumDeltas += fireworks[i].error - yMin;

      double[] result = new double[n]; // an amplitude for each firework
      double minAmplitude = MinAmplitude(epoch, maxEpochs, minX, maxX);
      for (int i = 0; i < n; ++i)
      {
        result[i] = A * (fireworks[i].error - yMin + 1.0E-10) / (sumDeltas + 1.0E-10);
        if (result[i] < minAmplitude)
          result[i] = minAmplitude;
      }
      return result;
    }

    private static double MinAmplitude(int epoch, int maxEpochs, double minX, double maxX)
    {
      // minimum amplitude for any firework at curr epoch 
      double Ainit = (maxX - minX) * 0.02;
      double Afinal = (maxX - minX) * 0.001;
      return Ainit - (Ainit - Afinal) * epoch / maxEpochs;
    }

    private static void AddGaussianSparks(Info[] fireworks, List<Info>[] sparksList, int dim, int mHat, int epoch, double minX, double maxX, double[] bestPosition, ref double bestError, Random rnd)
    {
      // generate mHat Gaussian sparks, add to sparksList, update bestPOsition[], bestError
      int n = fireworks.Length;
      for (int g = 0; g < mHat; ++g)
      {
        Info gSpark = new Info();
        gSpark.position = new double[dim];

        int i = rnd.Next(0, n); // pick a random firework
        for (int k = 0; k < dim; ++k) // spark position based on its parent firework
          gSpark.position[k] = fireworks[i].position[k];
        int z = (int)Math.Round(dim * rnd.NextDouble()); // number of random dimensions
        int[] dimensions = PickDimensions(dim, z, epoch); // pick random dimensions
        double u1 = rnd.NextDouble(); // make Gaussian displacement u = 0, sd = 1
        double u2 = rnd.NextDouble();
        double left = Math.Cos(2.0 * Math.PI * u1);
        double right = Math.Sqrt(-2.0 * Math.Log(u2));
        double e = left * right; // mean = 0, sd = 1
        for (int ii = 0; ii < dimensions.Length; ++ii) // each of the randomly selected dimensions
        {
          int k = dimensions[ii]; // convenience
          gSpark.position[k] = gSpark.position[k] + (bestPosition[k] - gSpark.position[k]) * e;
          if (gSpark.position[k] < minX || gSpark.position[k] > maxX) // bring out-of-range values back in
            gSpark.position[k] = (maxX - minX) * rnd.NextDouble() + minX;
        }
        gSpark.error = Error(gSpark.position);
        sparksList[i].Add(gSpark);

        if (gSpark.error < bestError)
        {
          bestError = gSpark.error;
          for (int k = 0; k < dim; ++k)
            bestPosition[k] = gSpark.position[k];
        }
      } // each Gaussian spark
    }

  } // Program

  public class Info
  {
    public double[] position;
    public double error;
  }

} // ns
