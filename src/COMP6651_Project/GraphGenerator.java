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
        int from, to;
        double x1, x2;
        double y1, y2;

        Edge(int from, double x1, double y1, int to, double x2, double y2) {
            this.from = from;
            this.to = to;
            this.x1 = x1;
            this.y1 = y1;
            this.x2 = x2;
            this.y2 = y2;
        }
    }

    private List<Vertex> vertices;
    private List<Edge> edges;
    private Random rand = new Random();
    private int n;
    private double r;


    public GraphGenerator(int n, double r) {
        this.n = n;
        this.r = r;

        vertices = new ArrayList<>();
        edges = new ArrayList<>();
        generateGraph(n, r);
    }

    private void generateGraph(int n, double r) {
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
                    if (distance <= r && !edgeExists(i+1, j+1)) {
                        edges.add(new Edge(i+1, u.x, u.y, j+1, v.x, v.y));
                    }
                }
            }
        }
    }

    private boolean edgeExists(int u, int v) {
        return edges.stream().anyMatch(edge -> (edge.from == u && edge.to == v) || (edge.from == v && edge.to == u));
    }

    public void saveGraph(String filename, String extension) throws IOException {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(filename + "." + extension))) {
            // Write n, r, and upperCap as the first line in the file
            if(extension.equalsIgnoreCase("mtx")){
                bw.write(n + " " + n + " " + edges.size() + "\n");
            }
            if(extension.equalsIgnoreCase("edges")){
                bw.write(edges.size() + " " + n + " " + n + "\n");
            }

            for (Edge edge : edges) {
                String line = edge.from + " " + edge.x1 + " " + edge.y1 + " " +  edge.to + " " + edge.x2 + " " + edge.y2 + "\n";
                bw.write(line);
            }
        }
    }

    public static void main(String[] args) {
        int n = 500;
        double r = 0.057;

        GraphGenerator generator = new GraphGenerator(n, r);

        try {
            generator.saveGraph("n" + n + "r" + r, "edges");
            System.out.println("Graph saved");
        } catch (IOException e) {
            System.out.println("Error writing to file: " + e.getMessage());
        }
    }
}
