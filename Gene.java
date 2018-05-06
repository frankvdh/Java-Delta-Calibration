
import java.util.Arrays;
import java.util.Random;

/**
 * Genetic algorithm. Based on C source code received via UseNet in 1991 -- the
 * author didn't even give his name!!:
 *
 * This version will find the global (as opposed to a local) maximum of a
 * multivariate function. That function need not be continuous nor
 * differentiable. Both these aspects are plus points when compared with
 * gradient-oriented algorithms.
 *
 * @author frankv Frank van der Hulst. drifter.frank@gmail.com
 */
public class Gene {

    /* This is a program where we use a genetic algorithm to solve a nonlinear unconstrained maximisation problem. */
    static final int SIZE_POPULATION = 65000;
    /* each generation contains so many animals */
    static final double CONVERGED = 20;
    /* convergence criterion. */
    static double COSMIC_RAY_DENSITY = 0.00004;
    static Random rand = new Random();   // Seed random number generator

    // Two population arrays are used alternately, to avoid re-allocating memory
    static Animal[][] pops = new Animal[2][SIZE_POPULATION];
    static int[] numSex = new int[SIZE_POPULATION];
    static int current = 0;
    static String firmware = "Repetier";
    static String result;

    static Thread thread = new Thread(() -> {
        doGenetic(SIZE_POPULATION, CONVERGED);
    });

    static Delta getResults() {
        pops[current][0].delta.normaliseEndstopAdjustments(firmware);
        result = "SConverged!\n"
                + "Best = " + pops[current][0].toString();
        return pops[current][0].delta;

    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        getArgs(args);
        result = "EInterrupted";
        doGenetic(SIZE_POPULATION, CONVERGED);
        getResults();
        System.out.println(result.substring(1));
    }

    public static void doGenetic(int N, double converged) {
        /* setup populations at time 0 */
        firstGen();

        /* We now have one generation: N, x and y (incl probabilities).
	In the main loop of the GA, we're going to constantly
	manipulate two generations: current and new. */
        int time = 1;
        while (!Thread.currentThread().isInterrupted()) {
            printGeneration(time, pops[current]);
            double oldBest = pops[current][0].y;
            int newPop = (current + 1) & 1;
            breed(pops[current], pops[newPop]);
            time++;
            current = newPop;
            if (Math.abs(oldBest - pops[current][0].y) < converged) {
                System.out.println("" + Math.abs(oldBest - pops[current][0].y));
                break;
            }
        }
    }

    /* Algorithm from here on works for any scalar function eval(x). */
 /* This creates two vectors of animals... the first generation,
     and 'blank' animals in the next generation. All allocation of
     animal memory is done here */
    static void firstGen() {
        for (int i = 0; i < pops[0].length; i++) {
            pops[0][i] = new Animal();
            pops[1][i] = new Animal();
        }
        Arrays.sort(pops[0]);
    }

    /* Data structures: at time t, the population will be
characterised by three things:
	N, the number of animals
	x, a vector with N elements where x_i is the attribute of animal i
	y, a vector where y_i = f(x_i)
     */
    static double uniform(double min, double max) {
        return (min + rand.nextDouble() * (max - min));
    }

    /* The calcSexProbs operation is three passes over the table of
size N.  It does 2N double comparisons takes 4N (double) flops. */
    static int calcSexProbs(Animal[] pop) {
        double midPoint, min, max;
        long sum;
        Arrays.sort(pop);
        min = pop[pop.length - 1].y;
        max = pop[0].y;

        /* now replace by fitness as follows: let midPoint=(min+max)/2.
		<= midPoint --> 0
		max-->1.0
         */
        double range = max - min;
        midPoint = (max + min) / 2.0;
        sum = 0;
        int i;
        for (i = 0; pop[i].y > midPoint; i++) {
            // Reduce values to the range 0-10000
            numSex[i] = (int) Math.round((pop[i].y - min) / range * 50000);
            sum += numSex[i];
        }
        int avg_loc = i;

        for (; i < pop.length; i++) {
            numSex[i] = 0;
        }
        /* these are probabilities that add up to sum... reduce them so that they total to SIZE_POPULATION*2  */
        int totalSex = 0;
        sum /= (pop.length * 2);
        for (i = 0; i < avg_loc; i++) {
            numSex[i] = (int) ((numSex[i] + sum / 2) / sum);
            totalSex += numSex[i];
        }

        /* Now have calculated number of times each animal will have sex,
           these total to *about* SIZE_POPULATION*2. Increase/decrease to 
        exactly SIZE_POPULATION*2
         */
        while (totalSex < pop.length * 2) {
            numSex[rand.nextInt(avg_loc)]++;
            totalSex++;
        }
        while (totalSex > pop.length * 2) {
            numSex[rand.nextInt(avg_loc)]--;
            totalSex--;
        }
        return avg_loc;
    }

    /* Create a new generation in newPop from the current generation . */
    static void breed(Animal[] current, Animal[] newPop) {
        int numBreeders = calcSexProbs(current);
        int p1 = 0;
        int i = 0;
        for (; i < newPop.length; p1++) {
            while (numSex[p1] > 0 && i < newPop.length) {
                // Choose a partner at random from those still available. Do not choose itself!
                // The last animal may be forced to breed with itself if all other animals
                // have already been paired
                int partner = (numBreeders == p1 + 1) ? p1 : p1 + 1 + rand.nextInt(numBreeders - p1 - 1);
                // If the randomly chosen partner isn't available, find the next worst available
                // This will always find a partner, since the random number <= last available
                while (numSex[partner] == 0 && partner < numBreeders) {
                    partner++;
                }
                numSex[partner]--;
                numSex[p1]--;
                // If this was the last available, reduce the number of available animals
                while (numBreeders > 0 && numSex[numBreeders - 1] == 0) {
                    numBreeders--;
                }
                current[p1].sex(current[partner], newPop[i++]);
            }
        }
    }

    static void printGeneration(int time, Animal[] p) {
        System.out.printf("\nGeneration %d, (%7.3f - %7.3f)\n", time, p[p.length - 1].y, p[0].y);
        int N = (p.length > 9) ? 9 : p.length;

        for (int i = 0; i < N; i++) {
            System.out.println(p[i].toString());
        }
    }

    /**
     * Command line syntax:
     *
     * GENE
     * [/P=pop_size][/R=rays][/B=bottom][/T=top][/C=convergence][/V=verbose]
     *
     * where pop_size = max. population size default = 20 rays = cosmic ray
     * density 0.4 bottom, top = limits for first generation -30.0, +30.0
     * verbose = 0 for no verbose output 0 convergence = convergence criterion
     * 0.001
     */
    static void getArgs(String[] args) {
        for (int i = 0; i < args.length; i++) {
            args[i] = args[i].toUpperCase();
            if (args[i].charAt(0) != '/' && args[i].charAt(0) != '-') {
                System.out.printf("Invalid command line switch: %s\n(Must start with '-' or '/')\n", args[i]);
                System.out.printf("Use /? to get help\n");
                System.exit(1);
            }
            if (args[i].charAt(1) == '?' || args[i].charAt(1) == 'H') {
                System.out.printf("Command syntax: %s [/P=pop_size][/R=rays][/B=bottom][/T=top][/C=convergence]\n\n", args[0]);
                System.out.printf("pop_size    = max. population size\n");
                System.out.printf("rays        = cosmic ray density\n");
                System.out.printf("convergence = convergence criterion\n\n");
                System.exit(0);
            }
            if (args[i].charAt(2) != '=') {
                System.out.printf("Invalid command line switch: %s\n(Must be /%c=VALUE)\n", args[i], args[i].charAt(1));
                System.out.printf("Use /? to get help\n");
                System.exit(1);
            }
            switch (args[i].charAt(1)) {
                case 'P':
                    //                   SIZE_POPULATION = Integer.parseInt(args[i].substring(3));
                    break;
                case 'C':
//                    CONVERGED = Double.parseDouble(args[i].substring(3));
                    break;
                default:
                    System.out.printf("Invalid command line switch: %s\n(Must be /P, /R, or /C)\n", args[i]);
                    System.out.printf("Use /? to get help\n");
                    System.exit(1);
            }
        }
    }
}
