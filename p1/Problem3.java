package project1;

public class Problem3 {
  private static final int NUM_PROCESSORS = 5;
  private static volatile int current = 1;
  private Object lock = new Object();
  ProcessorThread[] threads;

  public static void main(String[] args) {
    Problem3 problem = new Problem3();
    problem.run();
  }

  private void run() {
    createThreads();
    runThreads();
    waitAllThreads();
  }

  private void createThreads() {
    threads = new ProcessorThread[NUM_PROCESSORS];
    for (int i = 0; i < NUM_PROCESSORS; i++) {
      threads[i] = new ProcessorThread(i + 1);
    }
  }

  private void runThreads() {
    for (int i = NUM_PROCESSORS - 1; i >= 0; i--) {
      threads[i].start();
    }
  }

  private void waitAllThreads() {
    for (int i = 0; i < NUM_PROCESSORS; i++) {
      while (threads[i].isAlive()) {
        try {
          threads[i].join();
        } catch (InterruptedException e) {
        }
      }
    }
    System.out.println("END");
  }

  private class ProcessorThread extends Thread {
    int id;

    public ProcessorThread(int id) {
      this.id = id;
    }

    @Override
    public void run() {
      while (true) {
        if (this.id == current) {
          System.out.println(this.id);
          synchronized (lock) {
            current++;
          }
          return;
        } else {
          try {
            Thread.sleep(500);
          } catch (InterruptedException e) {
          }
        }
      }
    }
  }
}
