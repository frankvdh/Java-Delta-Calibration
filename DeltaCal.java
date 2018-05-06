
import mikera.vectorz.Vector3;

/**
 *
 * @author frank
 */
public class DeltaCal {
// Delta calibration script

    boolean debug = true;

    /**
     * @param args the command line arguments
     * @throws java.lang.Exception
     */
    public static void main(String args[]) throws Exception {
        Delta p = new Delta(235, 143, 300, 0, 0, 0, 210, 330, 90);
        double[] carriages = p.inverseKinematics(new Vector3(0, 0, 0));
        System.out.println(printVector("Carriages", carriages, 2));
        Vector3 pos = p.forwardKinematics(carriages);
        System.out.println(pos.toString());
    }

    DeltaCal() {
    }

    static String printVector(String label, double[] v, int dp) {
        String rslt = label + ": {";
        for (int i = 0; i < v.length; ++i) {
            rslt += String.format("%." + dp + 'f', v[i]);
            if (i + 1 != v.length) {
                rslt += ", ";
            }
        }
        rslt += "}";
        return rslt;
    }

    static String printTable(String label, double[][] v, int dp) {
        String rslt = label + ":\n";
        for (int i = 0; i < v.length; ++i) {
            for (int j = 0; j < v[i].length; j++) {
                rslt += String.format("%." + dp + 'f', v[i][j]);
                rslt += (j == v[i].length - 1) ? "\n" : ", ";
            }
        }
        return rslt;
    }

    String dc42Calculations(Delta delta, int numPoints, int factors, DoubleField[][] probes, String firmware, boolean normalise) {
        // Transform the probing points to motor endpoints and store them in a matrix,
        // so that we can do multiple iterations using the same data
        double[] zProbeValues = new double[numPoints];
        double[][] probeMotorPositions = new double[numPoints][3];
        double[] expectedResiduals = new double[numPoints];
        double initialSumOfSquares = 0.0;
        for (int i = 0; i < numPoints; ++i) {
            double[] machinePos = new double[3];
            machinePos[0] = probes[i][0].get();
            machinePos[1] = probes[i][1].get();
            machinePos[2] = 0.0;
            zProbeValues[i] = probes[i][2].get();
            initialSumOfSquares += Delta.square(zProbeValues[i]);

            probeMotorPositions[i][0] = delta.transform(machinePos, 0);
            probeMotorPositions[i][1] = delta.transform(machinePos, 1);
            probeMotorPositions[i][2] = delta.transform(machinePos, 2);
        }

        double expectedRmsError;
        try {
            expectedRmsError = doDeltaCalibration(numPoints, factors, firmware, delta, normalise, zProbeValues, probeMotorPositions, expectedResiduals);
            for (int i = 0; i < numPoints; i++) {
                probes[i][3].set(expectedResiduals[i], 4);
            }
        } catch (Exception err) {
            return ("E" + err.getMessage());
        }
        return String.format("SCalibrated %d factors using %d points\n"
                + "Deviation before %.2f, after %.2f", factors, numPoints, Math.sqrt(initialSumOfSquares / numPoints), expectedRmsError);
    }

    double doDeltaCalibration(int numPoints, int numFactors, String firmware,
            Delta params, boolean normalise, double[] zProbeValues,
            double[][] probeMotorPositions, double[] expectedResiduals) throws Exception {
        System.out.println(printTable("Motor positions", probeMotorPositions, 2));
        double[] corrections = new double[numPoints];
        // Do 1 or more Newton-Raphson iterations
        // Two is slightly better than one, but three doesn't improve things.
        // Alternatively, we could stop when the expected RMS error is only slightly worse than the RMS of the residuals.
        double expectedRmsError = Double.NaN;
        for (int iteration = 0; iteration < 2; iteration++) {
            // Build a Nx7 matrix of derivatives with respect to xa, xb, yc, za, zb, zc, diagonal.
            double[][] derivative = new double[numPoints][numFactors];
            for (int i = 0; i < numPoints; ++i) {
                for (int j = 0; j < numFactors; ++j) {
                    try {
                    derivative[i][j] = params.computeDerivative(j, probeMotorPositions[i][0], probeMotorPositions[i][1], probeMotorPositions[i][2]);
                    } catch(Exception ex) {
                        System.out.println("Derivative failed at point " + i + ", factor " + j);
                derivative[i][j] = Double.NaN;
                }
                }
            }

            System.out.println(printTable("Derivative matrix", derivative, 4));

            // Now build the normal equations for least squares fitting
            Matrix normalMatrix = new Matrix(numFactors, numPoints, derivative, zProbeValues, corrections);

            System.out.println(normalMatrix.print("Normal matrix:"));

            double[] solution = normalMatrix.gaussJordan(numFactors);
            for (int i = 0; i < numFactors; ++i) {
                if (Double.isNaN(solution[i])) {
                    throw new Exception("Unable to calculate corrections. Please make sure the bed probe points are all distinct.");
                }
            }
            System.out.println(normalMatrix.print("Solved matrix:"));

            if (debug) {
                System.out.println(printVector("Solution", solution, 4));

                // Calculate and display the residuals
                double[] residuals = new double[numPoints];
                for (int i = 0; i < numPoints; ++i) {
                    double r = zProbeValues[i];
                    for (int j = 0; j < numFactors; ++j) {
                        r += solution[j] * derivative[i][j];
                    }
                    residuals[i] = r;
                }
                System.out.println(printVector("Residuals", residuals, 4));
            }

            params.adjust(numFactors, solution, normalise, firmware);

            // Calculate the expected probe heights using the new parameters
            double sumOfSquares = 0.0;
            for (int i = 0; i < numPoints; ++i) {
                for (int axis = 0; axis < 3; ++axis) {
                    probeMotorPositions[i][axis] += solution[axis];
                }
                try {
                    double newZ = params.inverseTransform(probeMotorPositions[i]);
                    corrections[i] = newZ;
                    expectedResiduals[i] = zProbeValues[i] + newZ;
                    sumOfSquares += Delta.square(expectedResiduals[i]);
                } catch (Exception ex) {
                    corrections[i] = 0;
                    expectedResiduals[i] = Double.NaN;
                }
            }
            expectedRmsError = Math.sqrt(sumOfSquares / numPoints);
            System.out.println(printVector("Expected probe error", expectedResiduals, 4));
        }
        return expectedRmsError;
    }

    public static double square(double x) {
        return x * x;
    }
}
