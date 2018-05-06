
import java.util.logging.Level;
import java.util.logging.Logger;
import mikera.vectorz.Vector2;
import mikera.vectorz.Vector3;

/**
 *
 * @author frank
 */
public class Delta {

    double diagonal, radius, homedHeight, homedCarriageHeight, endstop[], towerAngle[];

    private double D2, Q, Q2;
    private Vector2 bc, ca, ab;
    private final double[] coreF = new double[3];
    private Vector2[] tower = new Vector2[3];

    final double Hez = 0;
    final double effectorRadius = 25;

    static final double DEGREES_TO_RADIANS = Math.PI / 180.0;

    Delta(double diagonal, double radius, double height,
            double xstop, double ystop, double zstop,
            double xTowerAngle, double yTowerAngle, double zTowerAngle) {
        this.diagonal = diagonal;
        this.radius = radius;
        this.homedHeight = height;
        endstop = new double[3];
        this.endstop[1] = xstop;
        this.endstop[2] = ystop;
        this.endstop[0] = zstop;
        towerAngle = new double[3];
        towerAngle[0] = xTowerAngle;
        towerAngle[1] = yTowerAngle;
        towerAngle[2] = zTowerAngle;
        try {
            recalc();
        } catch (Exception ex) {
            Logger.getLogger(Delta.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    final void recalc() throws Exception {
        D2 = square(diagonal);
        tower[0] = new Vector2(radius * Math.cos(towerAngle[0] * DEGREES_TO_RADIANS), radius * Math.sin(towerAngle[0] * DEGREES_TO_RADIANS));
        tower[1] = new Vector2(radius * Math.cos(towerAngle[1] * DEGREES_TO_RADIANS), radius * Math.sin(towerAngle[1] * DEGREES_TO_RADIANS));
        tower[2] = new Vector2(radius * Math.sin(towerAngle[2] * DEGREES_TO_RADIANS), radius * Math.cos(towerAngle[2] * DEGREES_TO_RADIANS));
        bc = new Vector2(tower[2].x - tower[1].x, tower[2].y - tower[1].y);
        ca = new Vector2(tower[0].x - tower[2].x, tower[0].y - tower[2].y);
        ab = new Vector2(tower[1].x - tower[0].x, tower[1].y - tower[1].y);
        for (int i = 0; i < 3; i++) {
            coreF[i] = square(tower[i].x) + square(tower[i].y);
        Q = 2 * (ca.x * ab.y - ab.x * ca.y);
        Q2 = square(Q);
        }

        // Calculate the base carriage height when the printer is homed.
        double tempHeight = diagonal;		// any sensible height will do here, probably even zero
        homedCarriageHeight = homedHeight; // + tempHeight - inverseTransform(tempHeight, tempHeight, tempHeight);
    }

    double computeDerivative(int deriv, double ha, double hb, double hc) throws Exception {
        double perturb = 0.2;			// perturbation amount in mm or degrees
        Delta hiParams = new Delta(diagonal, radius, homedHeight, endstop[1], endstop[2], endstop[0], towerAngle[1], towerAngle[2], towerAngle[0]);
        Delta loParams = new Delta(diagonal, radius, homedHeight, endstop[1], endstop[2], endstop[0], towerAngle[1], towerAngle[2], towerAngle[0]);
        switch (deriv) {
            case 0:
            case 1:
            case 2:
                break;

            case 3:
                hiParams.radius += perturb;
                loParams.radius -= perturb;
                break;

            case 4:
                hiParams.towerAngle[1] += perturb;
                loParams.towerAngle[1] -= perturb;
                break;

            case 5:
                hiParams.towerAngle[2] += perturb;
                loParams.towerAngle[2] -= perturb;
                break;

            case 6:
                hiParams.diagonal += perturb;
                loParams.diagonal -= perturb;
                break;
        }

        hiParams.recalc();
        loParams.recalc();

        double[] Hhi = {(deriv == 0) ? ha + perturb : ha, (deriv == 1) ? hb + perturb : hb, (deriv == 2) ? hc + perturb : hc};
        double[] Hlo = {(deriv == 0) ? ha - perturb : ha, (deriv == 1) ? hb - perturb : hb, (deriv == 2) ? hc - perturb : hc};
        double zHi = hiParams.inverseTransform(Hhi);
        double zLo = loParams.inverseTransform(Hlo);

        return (zHi - zLo) / (2 * perturb);
    }

    // Make the average of the endstop adjustments zero, or make all endstop 
    //corrections negative, without changing the individual homed carriage heights
    void normaliseEndstopAdjustments(String firmware) {
        double eav = (firmware.startsWith("Marlin") || firmware.equals("Repetier"))
                ? Math.min(endstop[0], Math.min(endstop[1], endstop[2]))
                : (endstop[0] + endstop[1] + endstop[2]) / 3.0;
        endstop[0] -= eav;
        endstop[1] -= eav;
        endstop[2] -= eav;
        homedHeight += eav;
        homedCarriageHeight += eav;		// no need for a full recalc, this is sufficient
    }

    // Inverse transform method, We only need the Z component of the result.
    double inverseTransform(double[] H) throws Exception {
        double[] F = {coreF[0] + square(H[0]), coreF[1] + square(H[1]), coreF[2] + square(H[2])};

        // Setup PQRSU such that x = -(S - uz)/P, y = (P - Rz)/Q
        double P = bc.x * F[0] + ca.x * F[1] + ab.x * F[2];
        double S = bc.y * F[0] + ca.y * F[1] + ab.y * F[2];

        double R = 2 * (bc.x * H[0] + ca.x * H[1] + ab.x * H[2]);
        double U = 2 * (bc.y * H[0] + ca.y * H[1] + ab.y * H[2]);

        double R2 = square(R), U2 = square(U);

        double A = U2 + R2 + Q2;
        double minusHalfB = S * U + P * R + H[0] * Q2 + tower[0].x * U * Q - tower[0].y * R * Q;
        double C = square(S + tower[0].x * Q) + square(P - tower[0].y * Q) + (square(H[0]) - D2) * Q2;

        double rslt = (minusHalfB - Math.sqrt(square(minusHalfB) - A * C)) / A;
        if (Double.isNaN(rslt)) {
            throw new Exception("At least one probe point is not reachable. Please correct your delta radius, diagonal rod length, or probe coordniates.");
        }
        return rslt;
    }

    // Perform 3, 4, 6 or 7-factor adjustment.
    // The input vector contains the following parameters in this order:
    //  X, Y and Z endstop adjustments
    //  If we are doing 4-factor adjustment, the next argument is the delta radius. Otherwise:
    //  X tower X position adjustment
    //  Y tower X position adjustment
    //  Z tower Y position adjustment
    //  Diagonal rod length adjustment
    void adjust(int numFactors, double[] v, boolean norm, String firmware) throws Exception {
        double oldCarriageHeightA = homedCarriageHeight + endstop[1];	// save for later

        // Update endstop adjustments
        endstop[0] += v[0];
        endstop[1] += v[1];
        endstop[2] += v[2];
        if (norm) {
            normaliseEndstopAdjustments(firmware);
        }

        if (numFactors >= 4) {
            radius += v[3];

            if (numFactors >= 6) {
                towerAngle[1] += v[4];
                towerAngle[2] += v[5];

                if (numFactors == 7) {
                    diagonal += v[6];
                }
            }
            recalc();
        }

        // Adjusting the diagonal and the tower positions affects the homed carriage height.
        // We need to adjust homedHeight to allow for this, to get the change that was requested in the endstop corrections.
        double heightError = homedCarriageHeight + endstop[1] - oldCarriageHeightA - v[0];
        homedHeight -= heightError;
        homedCarriageHeight -= heightError;
    }

    public double transform(double[] machinePos, int axis) {
        return machinePos[2] + Math.sqrt(D2 - square(machinePos[0] - tower[axis].x) - square(machinePos[1] - tower[axis].y));
    }

    public static double square(double x) {
        return x * x;
    }

    // From http://www.reprap.org/mediawiki/images/b/b5/Rostock_Delta_Kinematics_3.pdf
    /**
     *
     * @param dZ height of the carriage above the bed
     * @return [x,y,z] of intersection point; null if no intersection
     */
    public Vector3 forwardKinematics(double[] dZ) {
        // Get the carriage positions  
        Vector3[] carriagePos = {new Vector3(tower[0].x, tower[0].y, dZ[0]-endstop[0]),
            new Vector3(tower[1].x, tower[1].y, dZ[1]-endstop[0]),
            new Vector3(tower[2].x, tower[2].y, dZ[2]-endstop[0])};

        /* As discussed in https://en.wikipedia.org/wiki/Trilateration we are establishing a new coordinate system
    * in the plane of the three carriage points.
    *
    * This system will have the origin at carriagePos[0] and carriagePos[1] is on the x axis. 
    * carriagePos[2] is in the X-Y plane with a Z component of zero. 
    * We will define unit vectors in this coordinate system in our original coordinate system. Then, when we calculate
    * the x, y and z values, we can translate back into the original system by moving along those unit vectors
    * by the corresponding values.  We try to use the variable names from the Wikipedia article.
         */
        // Create a vector in old coords along x axis of new coord
        Vector3 ex = new Vector3(carriagePos[1]);
        ex.sub(carriagePos[0]);
        // Get the Magnitude of vector, & create unit vector.
        double d = ex.normalise();

        // Find vector from the origin of the new system to the third point
        Vector3 ey = new Vector3(carriagePos[2]);
        ey.sub(carriagePos[0]);
        // Use dot product to find the component of this vector on the X axis
        double i = ex.dotProduct(ey);
        // Create a vector along the x axis that represents the x component of p13
        Vector3 iex = ex.multiplyCopy(i);
        // Subtract the X component away from the original vector leaving only the Y component.
        ey.sub(iex);
        // Get the magnitude of Y component, and make a unit vector
        double j = ey.normalise();

        //The cross product of the unit x and y is the unit z
        Vector3 ez = new Vector3(ex);
        ez.crossProduct(ey);

        // Now we have the d, i, and j values defined in Wikipedia.
        // Plug them into the equations defined in Wikipedia for x, y and z
        //Xnew = (r1^2 - r2^2 – d^2)/2*d
        //     = (L^2 - L^2 - d^2)/2d = d/2; 
        double x = d / 2;
        //Ynew = (r1^2 - r3^2 + i^2 + j^2)/2*j - i * x/j
        //     = ((L^2 - L^2 + i^2 + j^2)/2 - i * x)/j
        //     = ((i^2 + j^2)/2 - i * x)/j
        double y = ((i * i + j * j) / 2 - i * x) / j;
        //Znew = sqrt(r1^2 - x^2 – y^2)
        //     = sqrt(L^2 - x^2 - y^2)
        double zSquared = D2 - square(x) - square(y);
        if (zSquared < 0) {
            // No real solution exists... spheres centred at the carriage positions don't intersect
            return null;
        }
        double z = Math.sqrt(zSquared);
        // Finally, start from the origin in the old coords and add vectors in the old coords that represent the
        // x, y, z to find the point in the old system
        Vector3 result = new Vector3(carriagePos[0]);
        result.add(ex.multiplyCopy(x));
        result.add(ey.multiplyCopy(y));
        result.add(ez.multiplyCopy(z));
        return result;
    }

    public double[] inverseKinematics(Vector3 target) {
        double[] result = new double[3];
        for (int i = 0; i < 3; i++) {
            Vector2 pivot = new Vector2(0, radius + effectorRadius).rotate2D(DEGREES_TO_RADIANS * towerAngle[i]);
            Vector2 effector = new Vector2(0, effectorRadius).rotate2D(DEGREES_TO_RADIANS * towerAngle[i]);
            Vector2 v = new Vector2(pivot.x - effector.x, pivot.y - effector.y);
            //  Now with our simplified model we have the following formulas.
            // (X - Avx)^2 + (Y - Avy)^2 = Ad^2 = D^2 - Acz^2
            // solve for Acz
            // Acz^2 = D^2 - (X - Avx)^2 - (Y - Avy)^2
            double cz = Math.sqrt(D2 - square(target.x - v.x) - square(target.y - v.y));
            // The real values we want is the distance of each carriage above the bed.
            // Subtract endstop so it gets included
            result[i] = target.z + cz + Hez;
        }
        return result;
    }

    @Override
    public String toString() {
        return String.format("%5.1f %5.1f %5.1f %5.2f %5.2f %5.2f %6.1f %6.1f %6.1f",
                diagonal, radius, homedHeight, endstop[0], endstop[1], endstop[2], towerAngle[0], towerAngle[1], towerAngle[2]);
    }
}
