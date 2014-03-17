package project2;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

public class Problem2 {
  public static void main(String args[]) {
    final String inputName = "adjmatrixB.txt";
    final String outputName = inputName.replace("adjmatrix", "eigenvalues");
    int[][] matrix = new int[32][32];
    try {
      FileInputStream fis = new FileInputStream(new File(inputName));
      InputStreamReader isr = new InputStreamReader(fis);
      BufferedReader br = new BufferedReader(isr);
      String line = null;
      int rowNumber = 0;
      while ((line = br.readLine()) != null) {
        int val;
        String[] cells = line.split(" ");
        for (int i = 0; i < 32; i++) {
          val = Integer.parseInt(cells[i]);
          matrix[rowNumber][i] = val;
        }
        rowNumber++;
      }
      br.close();
      double[][] dm = new double[32][32];
      for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 32; j++) {
          dm[i][j] = matrix[i][j];
        }
      }

      Matrix m = new Matrix(dm);
      EigenvalueDecomposition eig = m.eig();
      // Matrix ev = eig.getD();
      double[] img = m.eig().getImagEigenvalues();
      double[] real = m.eig().getRealEigenvalues();
      StringBuilder sb = new StringBuilder();
      for (double d : real) {
        sb.append(d);
        sb.append(System.lineSeparator());
      }

      FileOutputStream fos = new FileOutputStream(new File(outputName));
      OutputStreamWriter osr = new OutputStreamWriter(fos);
      BufferedWriter bw = new BufferedWriter(osr);

      bw.write(sb.toString());
      bw.close();

    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
