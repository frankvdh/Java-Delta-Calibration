
import mikera.vectorz.Vector3;

/**
 *
 * @author frank
 */
public class Animal implements Comparable<Animal> {

    Delta delta;
    double y;
    static final Delta NOMINAL = new Delta(235, 143, 300, 0, 0, 0, 210, 330, 90);
    static final double SPREAD = 0.1;

    static final Vector3[] MEASURED_POINTS = {
        new Vector3(0, 0, 0),
        new Vector3(0, 90, 0.72),
        new Vector3(-77.94, 45, -0.22),
        new Vector3(-77.94, -45, 0.3),
        new Vector3(0, -90, 0.35),
        new Vector3(77.94, -45, 1.9),
        new Vector3(77.94, 45, 3.2)
            };
    static double[][] motors = new double[MEASURED_POINTS.length][3];

    static {
        for (int i = 0; i < MEASURED_POINTS.length; i++) {
            Vector3 target = new Vector3(MEASURED_POINTS[i].x, MEASURED_POINTS[i].y, 0);
            motors[i] = NOMINAL.inverseKinematics(target);
        }

    }

    public Animal() {
        delta = new Delta(randomised(NOMINAL.diagonal, 230, 240), randomised(NOMINAL.radius, 135, 150), NOMINAL.homedHeight,
                randomised(NOMINAL.endstop[0], 0, 10), randomised(NOMINAL.endstop[1], 0, 10), randomised(NOMINAL.endstop[2], 0, 10),
                randomised(NOMINAL.towerAngle[0]), randomised(NOMINAL.towerAngle[1]), NOMINAL.towerAngle[2]);
        eval();
    }

    static final double randomised(double x) {
        return x == 0 ? (0.5 - Gene.rand.nextDouble()) * SPREAD : x + (0.5 - Gene.rand.nextDouble()) * x * SPREAD;
    }

    static final double randomised(double x, double min, double max) {
        return min +  Gene.rand.nextDouble() * (max-min);
    }

    /* Create a new child from the 2 parents. */
    void sex(Animal other, Animal child) {
        child.delta.diagonal = combine(delta.diagonal, other.delta.diagonal);
        child.delta.radius = combine(delta.radius, other.delta.radius);
        child.delta.homedHeight = delta.homedHeight;
        child.delta.endstop[0] = combine(delta.endstop[0], other.delta.endstop[0]);
        child.delta.endstop[1] = combine(delta.endstop[1], other.delta.endstop[1]);
        child.delta.endstop[2] = combine(delta.endstop[2], other.delta.endstop[2]);
        child.delta.towerAngle[0] = combine(delta.towerAngle[0], other.delta.towerAngle[0]);
        child.delta.towerAngle[1] = combine(delta.towerAngle[1], other.delta.towerAngle[1]);
        child.delta.towerAngle[2] = other.delta.towerAngle[2];
        child.eval();
    }

    static final double combine(double a, double b) {
        double rnd = (1 + (Gene.COSMIC_RAY_DENSITY * Gene.rand.nextGaussian()));
        if ((a+b) == 0)
            return rnd;
        return (a + b) / 2 * rnd;
    }

    static final double combineAngles(double a, double b) {
        return (a + b) / 2 * (1+ (Gene.COSMIC_RAY_DENSITY/1000 * Gene.rand.nextGaussian()));
    }

    /* Objective function */
 /* This function has two local maxima and one local minimum.
	f(0.6) = 3.8896
	f(1.85) = 0.8924
	f(3.35) = 5.8234	global maximum.
     */
//    final void eval() {
//        double x2, x3, x4;
//
//        x2 = x * x;
//        x3 = x2 * x;
//        x4 = x3 * x;
//        y = ((47 * x / 3) - (239 * x2 / 12) + (25 * x3 / 3) - (13 * x4 / 12));
//    }
    final void eval() {
        y = 0;
        for (int i = 0; i < motors.length; i++) {
            Vector3 pos = delta.forwardKinematics(motors[i]);
            y -= pos.distanceSquared(MEASURED_POINTS[i]);
        }
    }

    @Override
    public int compareTo(Animal o) {
        // Sort descending
        return Double.valueOf(o.y).compareTo(y);
    }

    @Override
    public String toString() {
        return delta.toString() + String.format(" y=%9.4f", y);
    }
}
