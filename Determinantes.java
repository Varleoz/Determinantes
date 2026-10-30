
import java.util.Locale;
import java.util.Random;
 
public class Determinantes {
 
    static final Random RNG = new Random(42); // semente fixa -> resultados reprodutíveis
 
    // =====================================================================
    // PARTE 1.1 - REGRA DE SARRUS (3x3), PASSO A PASSO
    // =====================================================================
    static double sarrus(double[][] m) {
        if (m.length != 3 || m[0].length != 3)
            throw new IllegalArgumentException("Sarrus só se aplica a matrizes 3x3.");
 
        // Diagonais principais (da esquerda para a direita)
        double d1 = m[0][0] * m[1][1] * m[2][2];
        double d2 = m[0][1] * m[1][2] * m[2][0];
        double d3 = m[0][2] * m[1][0] * m[2][1];
        // Diagonais secundárias (da direita para a esquerda)
        double s1 = m[0][2] * m[1][1] * m[2][0];
        double s2 = m[0][0] * m[1][2] * m[2][1];
        double s3 = m[0][1] * m[1][0] * m[2][2];
 
        double positivos = d1 + d2 + d3;
        double negativos = s1 + s2 + s3;
 
        System.out.println("  Passo 1 - diagonais principais:");
        System.out.printf("    a11*a22*a33 = %.2f%n    a12*a23*a31 = %.2f%n    a13*a21*a32 = %.2f%n", d1, d2, d3);
        System.out.printf("    Soma (positivos) = %.2f%n", positivos);
        System.out.println("  Passo 2 - diagonais secundárias:");
        System.out.printf("    a13*a22*a31 = %.2f%n    a11*a23*a32 = %.2f%n    a12*a21*a33 = %.2f%n", s1, s2, s3);
        System.out.printf("    Soma (negativos) = %.2f%n", negativos);
        System.out.printf("  Passo 3 - det = %.2f - (%.2f) = %.2f%n", positivos, negativos, positivos - negativos);
        return positivos - negativos;
    }
 
    // =====================================================================
    // DETERMINANTE VIA DECOMPOSIÇÃO LU COM PIVOTEAMENTO PARCIAL  O(n^3)
    // =====================================================================
    static double det(double[][] a) {
        int n = a.length;
        double[][] m = copia(a); // não altera a original
        double det = 1.0;
 
        for (int k = 0; k < n; k++) {
            // Pivoteamento parcial: maior |elemento| da coluna k -> estabilidade numérica
            int p = k;
            for (int i = k + 1; i < n; i++)
                if (Math.abs(m[i][k]) > Math.abs(m[p][k])) p = i;
 
            if (Math.abs(m[p][k]) < 1e-300) return 0.0; // matriz singular
 
            if (p != k) { // troca de linhas inverte o sinal do determinante
                double[] t = m[p]; m[p] = m[k]; m[k] = t;
                det = -det;
            }
            det *= m[k][k];
 
            // Eliminação: zera os elementos abaixo do pivô
            for (int i = k + 1; i < n; i++) {
                double f = m[i][k] / m[k][k];
                for (int j = k; j < n; j++) m[i][j] -= f * m[k][j];
            }
        }
        return det;
    }
 
    // =====================================================================
    // DETERMINANTE POR EXPANSÃO DE LAPLACE (RECURSIVO)  O(n!)
    // Usado apenas para mostrar a limitação dos métodos recursivos.
    // =====================================================================
    static double detLaplace(double[][] a) {
        int n = a.length;
        if (n == 1) return a[0][0];
        if (n == 2) return a[0][0] * a[1][1] - a[0][1] * a[1][0];
 
        double soma = 0;
        for (int j = 0; j < n; j++) {
            // Menor complementar removendo linha 0 e coluna j
            double[][] menor = new double[n - 1][n - 1];
            for (int i = 1; i < n; i++) {
                int c = 0;
                for (int k = 0; k < n; k++) {
                    if (k == j) continue;
                    menor[i - 1][c++] = a[i][k];
                }
            }
            soma += ((j % 2 == 0) ? 1 : -1) * a[0][j] * detLaplace(menor);
        }
        return soma;
    }
 
    // =====================================================================
    // INVERSA POR GAUSS-JORDAN COM PIVOTEAMENTO PARCIAL  O(n^3)
    // Retorna null se a matriz for singular.
    // =====================================================================
    static double[][] inversa(double[][] a) {
        int n = a.length;
        double[][] m = copia(a);
        double[][] inv = identidade(n);
 
        for (int k = 0; k < n; k++) {
            int p = k;
            for (int i = k + 1; i < n; i++)
                if (Math.abs(m[i][k]) > Math.abs(m[p][k])) p = i;
            if (Math.abs(m[p][k]) < 1e-12) return null; // singular
 
            double[] t = m[p]; m[p] = m[k]; m[k] = t;
            t = inv[p]; inv[p] = inv[k]; inv[k] = t;
 
            // Normaliza a linha do pivô
            double piv = m[k][k];
            for (int j = 0; j < n; j++) { m[k][j] /= piv; inv[k][j] /= piv; }
 
            // Zera a coluna k nas demais linhas
            for (int i = 0; i < n; i++) {
                if (i == k) continue;
                double f = m[i][k];
                if (f == 0) continue;
                for (int j = 0; j < n; j++) {
                    m[i][j] -= f * m[k][j];
                    inv[i][j] -= f * inv[k][j];
                }
            }
        }
        return inv;
    }
 
    // =====================================================================
    // UTILITÁRIOS DE MATRIZES
    // =====================================================================
    static double[][] copia(double[][] a) {
        double[][] c = new double[a.length][];
        for (int i = 0; i < a.length; i++) c[i] = a[i].clone();
        return c;
    }
 
    static double[][] identidade(int n) {
        double[][] I = new double[n][n];
        for (int i = 0; i < n; i++) I[i][i] = 1.0;
        return I;
    }
 
    static double[][] aleatoriaInteira(int n, int min, int max) {
        double[][] a = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                a[i][j] = min + RNG.nextInt(max - min + 1);
        return a;
    }
 
    static double[][] aleatoriaReal(int n) {
        double[][] a = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                a[i][j] = RNG.nextDouble() * 2 - 1; // [-1, 1)
        return a;
    }
 
    static double[][] multiplica(double[][] a, double[][] b) {
        int n = a.length;
        double[][] c = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++) {
                double aik = a[i][k];
                for (int j = 0; j < n; j++) c[i][j] += aik * b[k][j];
            }
        return c;
    }
 
    static double[][] transposta(double[][] a) {
        int n = a.length;
        double[][] t = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) t[j][i] = a[i][j];
        return t;
    }
 
    /** Maior |A*A^-1 - I|: mede a qualidade numérica da inversa. */
    static double erroInversa(double[][] a, double[][] inv) {
        double[][] p = multiplica(a, inv);
        double max = 0;
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a.length; j++)
                max = Math.max(max, Math.abs(p[i][j] - (i == j ? 1 : 0)));
        return max;
    }
 
    static void imprime(String titulo, double[][] a) {
        System.out.println(titulo);
        for (double[] linha : a) {
            StringBuilder sb = new StringBuilder("   [");
            for (int j = 0; j < linha.length; j++)
                sb.append(String.format(Locale.US, "%10.4f", linha[j]));
            System.out.println(sb.append(" ]"));
        }
    }
 
    static boolean iguais(double x, double y) {
        // tolerância relativa + absoluta
        return Math.abs(x - y) <= 1e-9 * Math.max(1.0, Math.max(Math.abs(x), Math.abs(y)));
    }
 
    // =====================================================================
    // PARTE 1 - EXECUÇÃO
    // =====================================================================
    static void parte1() {
        System.out.println("==================== PARTE 1 ====================");
 
        // ---- 1.1 Sarrus ----
        System.out.println("\n--- 1.1 Regra de Sarrus (3x3) ---");
        double[][] a3 = {{2, -1, 3}, {0, 4, 1}, {5, 2, -2}};
        imprime("Matriz A (3x3):", a3);
        double dS = sarrus(a3);
        System.out.printf(Locale.US, "  Conferência com LU: %.6f%n", det(a3));
        System.out.println("  Sarrus == LU ? " + iguais(dS, det(a3)));
 
        // ---- 1.2 / 1.3 / 1.4 Matrizes de teste ----
        System.out.println("\n--- 1.2 a 1.4 Determinante e inversa para n = 2, 3, 10, 50, 100 ---");
        int[] ordens = {2, 3, 10, 50, 100};
 
        for (int n : ordens) {
            System.out.println("\n>>> n = " + n);
            double[][] a;
            Double esperado = null;
 
            if (n <= 10) {
                // Ordens pequenas: inteiros aleatórios entre -9 e 9
                a = aleatoriaInteira(n, -9, 9);
            } else {
                // Ordens grandes: matriz diagonal com det conhecido analiticamente.
                // d_i = 1 para i par, 2 para i ímpar  =>  det = 2^(n/2)
                a = new double[n][n];
                int quantosDois = 0;
                for (int i = 0; i < n; i++) {
                    a[i][i] = (i % 2 == 0) ? 1.0 : 2.0;
                    if (i % 2 != 0) quantosDois++;
                }
                esperado = Math.pow(2, quantosDois);
            }
 
            if (n <= 3) imprime("Matriz A:", a);
 
            double d = det(a);
            System.out.printf(Locale.US, "det(A) [LU]            = %.6e%n", d);
            if (n == 3) System.out.printf(Locale.US, "det(A) [Sarrus]        = %.6e%n", sarrus3Silencioso(a));
            if (n <= 10) System.out.printf(Locale.US, "det(A) [Laplace]       = %.6e%n", detLaplace(a));
            if (esperado != null) {
                System.out.printf(Locale.US, "det(A) [analítico]     = %.6e  (erro relativo: %.2e)%n",
                        esperado, Math.abs(d - esperado) / esperado);
            }
 
            double[][] inv = inversa(a);
            if (inv == null) {
                System.out.println("Matriz singular: não possui inversa.");
            } else {
                if (n <= 3) imprime("Inversa A^-1:", inv); // só imprime inversas pequenas
                else System.out.println("(inversa omitida na tela por ser grande)");
                System.out.printf(Locale.US, "max|A*A^-1 - I|     = %.3e%n", erroInversa(a, inv));
            }
        }
 
        // Matriz identidade grande (det = 1) como segundo caso estruturado
        double[][] I = identidade(100);
        System.out.printf(Locale.US, "%n>>> Identidade 100x100: det = %.1f (esperado 1.0)%n", det(I));
    }
 
    /** Sarrus sem imprimir passos (apenas para conferência). */
    static double sarrus3Silencioso(double[][] m) {
        return m[0][0]*m[1][1]*m[2][2] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1]
             - m[0][2]*m[1][1]*m[2][0] - m[0][0]*m[1][2]*m[2][1] - m[0][1]*m[1][0]*m[2][2];
    }
 
    // =====================================================================
    // PARTE 2 - PROPRIEDADES (n = 5)
    // =====================================================================
    static void teste(String nome, double esq, double dir) {
        System.out.printf(Locale.US, "  [%s] %-52s %.6e vs %.6e%n",
                iguais(esq, dir) ? "OK " : "ERRO", nome, esq, dir);
    }
 
    static void parte2() {
        System.out.println("\n==================== PARTE 2 ====================");
        System.out.println("Propriedades verificadas com matrizes 5x5 (inteiros em [-5,5])\n");
        int n = 5;
        double[][] A = aleatoriaInteira(n, -5, 5);
        double[][] B = aleatoriaInteira(n, -5, 5);
        double dA = det(A), dB = det(B);
        System.out.printf(Locale.US, "det(A) = %.4f   det(B) = %.4f%n%n", dA, dB);
 
        // P1: det(A^T) = det(A)
        teste("P1: det(A^T) = det(A)", det(transposta(A)), dA);
 
        // P2: det(A*B) = det(A)*det(B)  (Teorema de Binet)
        teste("P2: det(AB) = det(A)*det(B)", det(multiplica(A, B)), dA * dB);
 
        // P3: det(I) = 1
        teste("P3: det(I) = 1", det(identidade(n)), 1.0);
 
        // P4: trocar duas linhas troca o sinal
        double[][] T = copia(A);
        double[] tmp = T[0]; T[0] = T[3]; T[3] = tmp;
        teste("P4: troca de duas linhas => det muda de sinal", det(T), -dA);
 
        // P5: multiplicar uma linha por k multiplica o det por k
        double k = 3.0;
        double[][] M = copia(A);
        for (int j = 0; j < n; j++) M[2][j] *= k;
        teste("P5: linha * k => det * k", det(M), k * dA);
 
        // P6: det(kA) = k^n det(A)
        double[][] KA = copia(A);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) KA[i][j] *= k;
        teste("P6: det(kA) = k^n * det(A)", det(KA), Math.pow(k, n) * dA);
 
        // P7: duas linhas iguais => det = 0
        double[][] L = copia(A);
        L[4] = L[1].clone();
        teste("P7: duas linhas iguais => det = 0", det(L), 0.0);
 
        // P8: linha nula => det = 0
        double[][] Z = copia(A);
        Z[2] = new double[n];
        teste("P8: linha nula => det = 0", det(Z), 0.0);
 
        // P9: somar múltiplo de uma linha a outra não altera o det (Jacobi)
        double[][] J = copia(A);
        for (int j = 0; j < n; j++) J[1][j] += 2.5 * J[3][j];
        teste("P9: L1 <- L1 + c*L3 => det inalterado", det(J), dA);
 
        // P10: det(A^-1) = 1/det(A)
        double[][] inv = inversa(A);
        if (inv != null) teste("P10: det(A^-1) = 1/det(A)", det(inv), 1.0 / dA);
        else System.out.println("  [--] P10: A singular, sem inversa.");
 
        // P11: matriz triangular superior => det = produto da diagonal
        double[][] U = copia(A);
        double prod = 1;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) U[i][j] = 0;
            if (U[i][i] == 0) U[i][i] = i + 1; // evita det trivialmente nulo
            prod *= U[i][i];
        }
        teste("P11: triangular => det = produto da diagonal", det(U), prod);
 
        // P12: linearidade em uma linha: det(..., u+v, ...) = det(..., u, ...) + det(..., v, ...)
        double[][] A1 = copia(A), A2 = copia(A), A3 = copia(A);
        for (int j = 0; j < n; j++) {
            A2[0][j] = RNG.nextInt(11) - 5;
            A3[0][j] = A1[0][j] + A2[0][j];
        }
        teste("P12: linearidade na linha 0", det(A3), det(A1) + det(A2));
    }
 
    // =====================================================================
    // PARTE 3 - DESEMPENHO
    // =====================================================================
    static void parte3() {
        System.out.println("\n==================== PARTE 3 ====================");
 
        // ---- 3.A: Laplace (recursivo) vs LU ----
        System.out.println("\n--- 3.A Expansão de Laplace O(n!) vs. LU O(n^3) ---");
        System.out.printf("%4s | %16s | %16s | %12s%n", "n", "Laplace (ms)", "LU (ms)", "Laplace/LU");
        System.out.println("-----+------------------+------------------+-------------");
        for (int n = 2; n <= 11; n++) {
            double[][] a = aleatoriaReal(n);
            long t0 = System.nanoTime();
            double dl = detLaplace(a);
            long t1 = System.nanoTime();
            double dlu = det(a);
            long t2 = System.nanoTime();
            double msL = (t1 - t0) / 1e6, msLU = Math.max((t2 - t1) / 1e6, 1e-6);
            System.out.printf(Locale.US, "%4d | %16.4f | %16.6f | %12.1f%s%n",
                    n, msL, msLU, msL / msLU, iguais(dl, dlu) ? "" : "  (valores divergem!)");
        }
 
        // ---- 3.B: LU - crescimento com n ----
        System.out.println("\n--- 3.B Tempo do determinante (LU) em função de n ---");
 
        // Aquecimento da JVM (JIT) para não distorcer as primeiras medições
        for (int i = 0; i < 30; i++) det(aleatoriaReal(200));
 
        int[] ordens = {10, 25, 50, 100, 200, 400, 800, 1600};
        double[] tempos = new double[ordens.length];
 
        System.out.printf("%6s | %14s | %18s | %14s%n", "n", "Tempo (ms)", "t(n)/t(n/2) ", "t/n^3 (ns)");
        System.out.println("-------+----------------+--------------------+---------------");
        for (int idx = 0; idx < ordens.length; idx++) {
            int n = ordens[idx];
            double[][] a = aleatoriaReal(n);
            int reps = (n <= 100) ? 200 : (n <= 400 ? 10 : (n <= 800 ? 3 : 1));
            double melhor = Double.MAX_VALUE; // menor tempo = menos ruído do SO/GC
            for (int r = 0; r < reps; r++) {
                long t0 = System.nanoTime();
                det(a);
                melhor = Math.min(melhor, (System.nanoTime() - t0) / 1e6);
            }
            tempos[idx] = melhor;
            String razao = "-";
            if (idx > 0 && ordens[idx] == 2 * ordens[idx - 1])
                razao = String.format(Locale.US, "%.2f", tempos[idx] / tempos[idx - 1]);
            System.out.printf(Locale.US, "%6d | %14.4f | %18s | %14.4f%n",
                    n, melhor, razao, melhor * 1e6 / ((double) n * n * n));
        }
 
        System.out.println("\nInterpretação:");
        System.out.println("  - Se T(n) = c*n^3, então dobrar n deve multiplicar o tempo por ~8 (coluna t(n)/t(n/2)).");
        System.out.println("  - A coluna t/n^3 deve ser aproximadamente constante (valor de c).");
        System.out.println("  - Para n pequeno há desvios (JIT, cache, overhead); para n grande,");
        System.out.println("    a hierarquia de cache pode elevar a razão acima de 8.");
 
        // ---- Saída CSV para gráfico ----
        System.out.println("\n--- CSV (para gráfico no Excel/Python/LibreOffice) ---");
        System.out.println("n,tempo_ms");
        for (int i = 0; i < ordens.length; i++)
            System.out.printf(Locale.US, "%d,%.6f%n", ordens[i], tempos[i]);
    }
 
    // =====================================================================
    public static void main(String[] args) {
        parte1();
        parte2();
        parte3();
    }
}
    