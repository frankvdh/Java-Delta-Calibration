
/**
 *
 * @author frank
 */
public class Matrix {

  private double[][] data;

  Matrix(int rows, int cols) {
    data = new double[rows][cols];
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        data[i][j] = 0.0;
      }
    }
  }

  public Matrix(int numFactors, int numPoints, double[][] derivative, double[] zBedProbePoints, double[] corrections) {
    // Now build the normal equations for least squares fitting
    this(numFactors, numFactors + 1);
    for (int i = 0; i < numFactors; ++i) {
      for (int j = 0; j < numFactors; ++j) {
        double temp = derivative[0][i] * derivative[0][j];
        for (int k = 1; k < numPoints; ++k) {
          temp += derivative[k][i] * derivative[k][j];
        }
        data[i][j] = temp;
      }
      double temp = derivative[0][i] * -(zBedProbePoints[0] + corrections[0]);
      for (int k = 1; k < numPoints; ++k) {
        temp += derivative[k][i] * -(zBedProbePoints[k] + corrections[k]);
      }
      data[i][numFactors] = temp;
    }
  }

  // Perform Gauus-Jordan elimination on a matrix with numRows rows and (numRows + 1) columns
  public double[] gaussJordan(int numRows) {
    double[] solution = new double[numRows];
    for (int i = 0; i < numRows; ++i) {
      // Swap the rows around for stable Gauss-Jordan elimination
      double vmax = Math.abs(data[i][i]);
      for (int j = i + 1; j < numRows; ++j) {
        double rmax = Math.abs(data[j][i]);
        if (rmax > vmax) {
          this.swapRows(i, j, numRows + 1);
          vmax = rmax;
        }
      }

      // Use row i to eliminate the ith element from previous and subsequent rows
      double v = data[i][i];
      for (int j = 0; j < i; ++j) {
        double factor = data[j][i] / v;
        data[j][i] = 0.0;
        for (int k = i + 1; k <= numRows; ++k) {
          data[j][k] -= data[i][k] * factor;
        }
      }

      for (int j = i + 1; j < numRows; ++j) {
        double factor = data[j][i] / v;
        data[j][i] = 0.0;
        for (int k = i + 1; k <= numRows; ++k) {
          data[j][k] -= data[i][k] * factor;
        }
      }
    }

    for (int i = 0; i < numRows; ++i) {
      solution[i] = data[i][numRows] / data[i][i];
    }
    return solution;
  }

  public void swapRows(int i, int j, int numCols) {
    if (i != j) {
      for (int k = 0; k < numCols; ++k) {
        double temp = data[i][k];
        data[i][k] = data[j][k];
        data[j][k] = temp;
      }
    }
  }

  public String print(String tag) {
    String rslt = tag + " \n";
    for (int i = 0; i < data.length; ++i) {
      double[] row = data[i];
      rslt += ' ';
      for (int j = 0; j < row.length; ++j) {
        rslt += String.format("%.4f", row[j]);
        if (j + 1 < row.length) {
          rslt += ", ";
        }
      }
      rslt += "\n";
    }
    return rslt;
  }
}
