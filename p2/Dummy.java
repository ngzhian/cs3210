package project2;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;

public class Dummy {
  public static void Dumb() throws IOException {
    /* read in adjlist and convert to adjmtx */
    final String inputName = "adjlistB.txt";
    final String outputName = inputName.replace("list", "matrix");
    int[][] matrix = new int[32][32];
    FileInputStream fis = new FileInputStream(new File(inputName));
    InputStreamReader isr = new InputStreamReader(fis);
    BufferedReader br = new BufferedReader(isr);
    String line = null;
    while ((line = br.readLine()) != null) {
      System.out.println("Read: " + line);
      // format is (x,y)
      String xcommay = line.substring(1, line.length() - 1);
      String[] coords = xcommay.split(",");
      int x = Integer.parseInt(coords[0]);
      int y = Integer.parseInt(coords[1]);
      // System.out.println("x: " + x + " y: " + y);
      matrix[x - 1][y - 1] = 1;
      matrix[y - 1][x - 1] = 1;
    }
    br.close();
    FileOutputStream fos = new FileOutputStream(new File(outputName));
    OutputStreamWriter osr = new OutputStreamWriter(fos);
    BufferedWriter bw = new BufferedWriter(osr);
    StringBuilder sb;
    for (int[] element : matrix) {
      sb = new StringBuilder();
      for (int val : element) {
        sb.append(val);
        sb.append(" ");
      }
      bw.write(sb.toString());
      bw.newLine();
    }
    bw.close();

    /* read in adj mtx */
    FileInputStream fis1 = new FileInputStream(new File(inputName));
    InputStreamReader isr1 = new InputStreamReader(fis1);
    BufferedReader br1 = new BufferedReader(isr1);
    String line1 = null;
    int rowNumber = 0;
    while ((line = br1.readLine()) != null) {
      int val;
      String[] cells = line.split(" ");
      for (int i = 0; i < 32; i++) {
        val = Integer.parseInt(cells[i]);
        matrix[rowNumber][i] = val;
      }
      rowNumber++;
    }
    br.close();

    /* helper to make adjlist for connected components */
    FileOutputStream fos2 = new FileOutputStream(new File("morelist.txt"));
    OutputStreamWriter osr2 = new OutputStreamWriter(fos2);
    BufferedWriter bw2 = new BufferedWriter(osr2);
    StringBuilder sb2;
    for (int i = 0; i < 8; i++) {
      int start = i * 4;
      sb2 = new StringBuilder();
      sb2.append(String.format("(%d,%d)", start, start + 1));
      sb2.append(System.lineSeparator());
      sb2.append(String.format("(%d,%d)", start + 1, start + 2));
      sb2.append(System.lineSeparator());
      sb2.append(String.format("(%d,%d)", start + 2, start + 3));
      sb2.append(System.lineSeparator());
      sb2.append(String.format("(%d,%d)", start, start + 3));
      sb2.append(System.lineSeparator());
      bw2.write(sb2.toString());
    }

    bw2.close();

  }

  public static int[][] fw(int mtx[][]) {
    int[][] ans = mtx.clone();
    int length = mtx[0].length;
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < length; j++) {
        if ((i != j) && (ans[i][j] == 0)) {
          ans[i][j] = 99999;
        }
      }
    }
    for (int k = 0; k < length; k++) {
      for (int i = 0; i < length; i++) {
        for (int j = 0; j < length; j++) {
          if (ans[i][j] > (ans[i][k] + ans[k][j])) {
            ans[i][j] = ans[i][k] + ans[k][j];
          }
        }
      }
    }
    return ans;
  }

  /* calculate hop distribution */
  public static int[] calcDistribution(int[][] matrix) {
    int[] distribution = new int[10];
    for (int i = 0; i < 32; i++) {
      for (int j = 0; j < 32; j++) {
        distribution[matrix[i][j]]++;
      }
    }
    return distribution;
  }

}