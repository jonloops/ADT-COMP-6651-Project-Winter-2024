package COMP6651_Project;

import java.io.*;
import java.util.*;

public class GraphGenerator {
    static class Vertex {
        double x, y;

        Vertex(double x, double y) {
            this.x = x;
            this.y = y;
        }
    }

    static class Edge {
        int from, to, capacity;

        Edge(int from, int to, int capacity) {
            this.from = from;
            this.to = to;
            this.capacity = capacity;
        }
    }

    private List<Vertex> vertices;
    private List<Edge> edges;
    private Random rand = new Random();
    private int n;
    private double r;
    private int upperCap;

    public GraphGenerator(int n, double r, int upperCap) {
        this.n = n;
        this.r = r;
        this.upperCap = upperCap;
        vertices = new ArrayList<>();
        edges = new ArrayList<>();
        generateGraph(n, r, upperCap);
    }

    private void generateGraph(int n, double r, int upperCap) {
        for (int i = 0; i < n; i++) {
            double x = rand.nextDouble();
            double y = rand.nextDouble();
            vertices.add(new Vertex(x, y));
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    Vertex u = vertices.get(i);
                    Vertex v = vertices.get(j);
                    double distance = Math.sqrt(Math.pow(u.x - v.x, 2) + Math.pow(u.y - v.y, 2));
                    if (distance <= r && !edgeExists(i, j)) {
                        int capacity = 1 + rand.nextInt(upperCap);
                        edges.add(new Edge(i, j, capacity));
                    }
                }
            }
        }
    }

    private boolean edgeExists(int u, int v) {
        return edges.stream().anyMatch(edge -> (edge.from == u && edge.to == v) || (edge.from == v && edge.to == u));
    }

    public void saveGraphToCSV(String filename) throws IOException {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(filename))) {
            // Write n, r, and upperCap as the first line in the file
            bw.write(n + "," + r + "," + upperCap + "\n");
            for (Edge edge : edges) {
                String line = edge.from + "," + edge.to + "," + edge.capacity + "\n";
                bw.write(line);
            }
        }
    }

    public static void main(String[] args) {
        int n = 7;
        double r = 0.4;
        int upperCap = 5;
        GraphGenerator generator = new GraphGenerator(n, r, upperCap);

        try {
            generator.saveGraphToCSV("graphTestManual.csv");
            System.out.println("Graph saved to graphTestManual.csv");
        } catch (IOException e) {
            System.out.println("Error writing to file: " + e.getMessage());
        }
    }
}
