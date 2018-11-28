/**
 *  Genetic Algorithm solution for a Floating point curve problem
 *  @Version 1.0
 *  @Author Nayra Adel
 **/

import javafx.util.Pair;

import java.util.*;

import static java.lang.Math.pow;

public class GA {

    private Scanner scanner;

    private final int population_size             = 1000;
    private final int maximum_generations         = 1000;
    private final double probability_of_crossover = 0.7;
    private final double probability_of_mutation  = 0.7;

    private final double upper_bound = 10;
    private final double lower_bound = -10;

    private int number_of_points;
    private int degree;

    private List<List<Double>> population;

    private List<Double> fitness;
    private List<Double> cumulative_Finesses;

    private List<Pair<Double, Double>> points;

    public static void main(String[] args) { new GA(); }

    private GA() {

        init();
        getUserInput();
        initPopulation();
        evaluatePopulationFitness();
        generations();
    }

    private void init() {

        scanner = new Scanner(System.in);

        number_of_points = 0;
        degree = 0;

        population           = new ArrayList<>();
        fitness              = new ArrayList<>();
        cumulative_Finesses  = new ArrayList<>();
        points               = new ArrayList<>();

        for (int i = 0; i < population_size; ++i) { fitness.add(0.0); }
    }

    private void getUserInput() {

        System.out.println("Enter the number of points: ");
        number_of_points = scanner.nextInt();

        System.out.println("Enter degree: ");
        degree = scanner.nextInt() + 1;

        System.out.println("Enter points: ");

        for (int j = 0; j < number_of_points; j++) {

            points.add(new Pair<>(scanner.nextDouble(), scanner.nextDouble()));
        }
    }

    private void generations() {

        for (int i = 0; i < maximum_generations; ++i) {

            nextGeneration();
            evaluatePopulationFitness();
        }
        printCoefficients(getMinFitness());
    }

    private void initPopulation() {

        for (int i = 0; i < population_size; ++i) population.add(generateChromosome());
    }

    // Generates a chromosome "Array of floating points from -10 to 10"
    private List generateChromosome() {

        List<Double> chromosome = new ArrayList<>();
        Double random;

        for (int j = 0; j < degree; ++j) {

            random = (Math.random() * (upper_bound - lower_bound)) + lower_bound;
            chromosome.add(random);
        }
        return chromosome;
    }

    // degrees 3 means => y-calc = a0*x^0 + a1*x^1 + a2*x^2 + a3*x^3
    // compare y-calc with ‘y’ elli kant m3 al ‘x’( Point(x, y))
    // Fitness Function
    private void evaluatePopulationFitness() {

        double error;
        for (int i = 0; i < population_size; ++i) {

            error = evaluateChromosomeFitness(population.get(i)) / number_of_points;
            fitness.set(i, error);
        }
        getMinFitness();
    }

    private double evaluateChromosomeFitness(List<Double> genes) {

        double y_calc = 0, gene, min_square_error = 0;

        for (int i = 0; i < number_of_points; ++i) {
            for (int j = 0; j < degree; j++) {

                gene = genes.get(j);
                y_calc += gene * pow(points.get(i).getKey(), j);
            }

            min_square_error += pow(y_calc - points.get(i).getValue(), 2);
        }
        return min_square_error;
    }

    // Selection
    private double calculate_Cumulative_Finesses() {

        double totalFitness = 0;
        for (int i = 0; i < population_size; ++i) {

            totalFitness += fitness.get(i);
            cumulative_Finesses.add(totalFitness);
        }
        return totalFitness;
    }

    private int rouletteWheelSelection() {

        double totalFitness = calculate_Cumulative_Finesses();
        double randomNum = Math.random() * (totalFitness + 1);
        double lowerLimit = 0;

        for (int i = 0; i < cumulative_Finesses.size(); ++i) {

            if (randomNum >= lowerLimit && randomNum < cumulative_Finesses.get(i)) return i;
            lowerLimit = cumulative_Finesses.get(i);
        }
        return 0;
    }

    // create a new generation's population then call the cross over
    private void nextGeneration() {

        int gene_1 = rouletteWheelSelection();
        int gene_2 = rouletteWheelSelection();

        while (gene_1 == gene_2) { gene_2 = rouletteWheelSelection(); }
        crossoverGenes(gene_1, gene_2);
    }

    // Crossover - Mutation - Replacement
    private void crossoverGenes(int parent1, int parent2) {

        List<Double> gene_1;
        List<Double> gene_2;

        double random_crossover = Math.random();
        int cross_point;

        if (random_crossover <= probability_of_crossover) { // r < PC

            cross_point = (int) ((Math.random() * (degree - 2) + 1) + 1);

            gene_1 = new ArrayList<>(population.get(parent1).subList(0, cross_point));
            gene_1.addAll(population.get(parent2).subList(cross_point, degree));

            gene_2 = new ArrayList<>(population.get(parent2).subList(0, cross_point));
            gene_2.addAll(population.get(parent1).subList(cross_point, degree));

            // Check mutation
            mutateGene(gene_1, parent1);
            mutateGene(gene_2, parent2);

            // Add new genes to population (replace his parent)
            population.set(parent1, gene_1);
            population.set(parent2, gene_2);
        } else {
            // Check mutation
            mutateGene(population.get(parent1), parent1);
            mutateGene(population.get(parent2), parent2);

            // Add the same parents after mutation to population
            population.set(parent1, population.get(parent1));
            population.set(parent2, population.get(parent2));
        }
    }

    private void mutateGene(List<Double> chromosome, int i) {

        double xi, y, dlxi, duxi, random, b, amount_of_mutation, r;
        for (int j = 0; j < degree; ++j) {

            xi = chromosome.get(j);

            dlxi = xi - lower_bound;
            duxi = upper_bound - xi;

            random = Math.random();
            double check_mutation = Math.random();

            if(check_mutation <= probability_of_mutation) {
                if (random <= 0.5) y = dlxi;
                else y = duxi;

                // di(xi, t) = y.(1 - r^(((1-t)/T)^b)
                r = Math.random();
                b = (Math.random() * ((5 - 0.5) + 1)) + 0.5; // b = dependency Factor
                amount_of_mutation = y * (1 - (pow(r, pow((1 - i) / maximum_generations, b))));

                chromosome.set(j, xi + amount_of_mutation);
            }
        }
    }

    private int getMinFitness() {

        return fitness.indexOf(Collections.min(fitness));
    }

    private void printCoefficients(int index){
        for(int i = 0 ; i < degree ; ++i)

            System.out.println("coff" + i + ": " + population.get(index).get(i));
    }
}