
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

// --- Konstanta dari Studi Kasus ---
const double L = 1.0;    // Induktansi (H)
const double R = 0.0;    // Resistansi (Ohm) - Sesuai kasus di paper
const double C = 0.25;   // Kapasitansi (F)
const double E0 = 1.0;   // Amplitudo Tegangan (V)
const double omega = 3.5;// Frekuensi Sudut (rad/s)

// --- Parameter Simulasi ---
const double t_max = 20.0; // Waktu simulasi maksimum (s)
const double h = 0.01;   // Ukuran langkah waktu (s)

// --- Definisi Sistem Persamaan Diferensial ---
// Model rangkaian RLC diubah menjadi sistem dua ODE orde 1:
// dy1/dt = y2
// dy2/dt = (1/L) * [E(t) - R*y2 - (1/C)*y1]
// y1 adalah muatan (q), y2 adalah arus (i)

/**
 * @brief Menghitung turunan dari muatan (dq/dt = i).
 * @param t Waktu saat ini (tidak digunakan di fungsi ini, tapi ada untuk konsistensi).
 * @param y1 Nilai muatan (q) saat ini.
 * @param y2 Nilai arus (i) saat ini.
 * @return Nilai dari dy1/dt.
 */
double f1(double t, double y1, double y2) {
    return y2;
}

/**
 * @brief Menghitung turunan dari arus (di/dt).
 * @param t Waktu saat ini.
 * @param y1 Nilai muatan (q) saat ini.
 * @param y2 Nilai arus (i) saat ini.
 * @return Nilai dari dy2/dt.
 */
double f2(double t, double y1, double y2) {
    double E_t = E0 * sin(omega * t); // Tegangan sumber sinusoidal
    return (1.0 / L) * (E_t - R * y2 - (1.0 / C) * y1);
}

/**
 * @brief Melakukan satu langkah perhitungan dengan metode Runge-Kutta orde 4.
 * Digunakan untuk inisialisasi metode Adams-Moulton.
 * @param t Waktu saat ini.
 * @param y1 Nilai y1 saat ini.
 * @param y2 Nilai y2 saat ini.
 * @param h Ukuran langkah.
 * @return Vektor berisi {y1_baru, y2_baru}.
 */
std::vector<double> rk4_step(double t, double y1, double y2, double h) {
    double k1_y1 = f1(t, y1, y2);
    double k1_y2 = f2(t, y1, y2);

    double k2_y1 = f1(t + 0.5 * h, y1 + 0.5 * h * k1_y1, y2 + 0.5 * h * k1_y2);
    double k2_y2 = f2(t + 0.5 * h, y1 + 0.5 * h * k1_y1, y2 + 0.5 * h * k1_y2);

    double k3_y1 = f1(t + 0.5 * h, y1 + 0.5 * h * k2_y1, y2 + 0.5 * h * k2_y2);
    double k3_y2 = f2(t + 0.5 * h, y1 + 0.5 * h * k2_y1, y2 + 0.5 * h * k2_y2);

    double k4_y1 = f1(t + h, y1 + h * k3_y1, y2 + h * k3_y2);
    double k4_y2 = f2(t + h, y1 + h * k3_y1, y2 + h * k3_y2);

    double y1_new = y1 + (h / 6.0) * (k1_y1 + 2.0 * k2_y1 + 2.0 * k3_y1 + k4_y1);
    double y2_new = y2 + (h / 6.0) * (k1_y2 + 2.0 * k2_y2 + 2.0 * k3_y2 + k4_y2);

    return {y1_new, y2_new};
}

/**
 * @brief Melakukan satu langkah dengan metode Adams-Bashforth/Moulton orde 4.
 * @param t Waktu di awal langkah (t_i).
 * @param y1 Nilai y1 di t_i.
 * @param y2 Nilai y2 di t_i.
 * @param h Ukuran langkah.
 * @param f_history Vektor berisi histori turunan dari 4 titik sebelumnya.
 * @return Vektor berisi nilai y1 dan y2 yang sudah dikoreksi di t_{i+1}.
 */
std::vector<double> adams_moulton_step(double t, double y1, double y2, double h, 
                                     const std::vector<std::vector<double>>& f_history) {
    
    // --- PREDICTOR (Adams-Bashforth 4th order) ---
    // y_p = y_i + (h/24) * [55*f_i - 59*f_{i-1} + 37*f_{i-2} - 9*f_{i-3}]
    // f_history[0] -> turunan di t_{i-3}, f_history[3] -> turunan di t_i
    double y1_p = y1 + (h / 24.0) * (55 * f_history[3][0] - 59 * f_history[2][0] + 37 * f_history[1][0] - 9 * f_history[0][0]);
    double y2_p = y2 + (h / 24.0) * (55 * f_history[3][1] - 59 * f_history[2][1] + 37 * f_history[1][1] - 9 * f_history[0][1]);

    // Hitung turunan di titik prediksi untuk digunakan di corrector
    double f1_p = f1(t + h, y1_p, y2_p);
    double f2_p = f2(t + h, y1_p, y2_p);

    // --- CORRECTOR (Adams-Moulton 4th order) ---
    // y_c = y_i + (h/24) * [9*f_{i+1} + 19*f_i - 5*f_{i-1} + f_{i-2}]
    // f_{i+1} adalah turunan di titik prediksi (f1_p, f2_p)
    double y1_c = y1 + (h / 24.0) * (9 * f1_p + 19 * f_history[3][0] - 5 * f_history[2][0] + f_history[1][0]);
    double y2_c = y2 + (h / 24.0) * (9 * f2_p + 19 * f_history[3][1] - 5 * f_history[2][1] + f_history[1][1]);

    return {y1_c, y2_c};
}


int main() {
    // Jumlah total iterasi
    int n_steps = static_cast<int>(t_max / h);

    // Vektor untuk menyimpan hasil simulasi (q dan i) sepanjang waktu
    std::vector<double> time(n_steps + 1);
    std::vector<double> q_charge(n_steps + 1); // y1
    std::vector<double> i_current(n_steps + 1); // y2

    // Mengatur kondisi awal: q(0) = 0, i(0) = 0
    time[0] = 0.0;
    q_charge[0] = 0.0;
    i_current[0] = 0.0;

    // Vektor untuk menyimpan histori turunan [f1, f2] yang dibutuhkan metode Adams
    std::vector<std::vector<double>> f_history(4, std::vector<double>(2));

    // --- FASE INISIALISASI DENGAN RK4 ---
    // Metode Adams butuh 4 titik awal (t0, t1, t2, t3) untuk memulai.
    // Kita gunakan RK4 untuk menghasilkan 3 titik berikutnya setelah t0.
    
    // Hitung turunan di titik awal t=0
    f_history[0] = {f1(time[0], q_charge[0], i_current[0]), f2(time[0], q_charge[0], i_current[0])};

    for (int i = 0; i < 3; ++i) {
        std::vector<double> next_y = rk4_step(time[i], q_charge[i], i_current[i], h);
        time[i+1] = time[i] + h;
        q_charge[i+1] = next_y[0];
        i_current[i+1] = next_y[1];
        
        // Hitung dan simpan histori turunan untuk digunakan metode Adams
        f_history[i+1] = {f1(time[i+1], q_charge[i+1], i_current[i+1]), f2(time[i+1], q_charge[i+1], i_current[i+1])};
    }

    // --- FASE UTAMA DENGAN ADAMS-MOULTON ---
    // Loop utama dimulai dari titik ke-3 karena 3 titik pertama sudah dihitung oleh RK4
    for (int i = 3; i < n_steps; ++i) {
        time[i+1] = time[i] + h;
        
        // Lakukan satu langkah perhitungan dengan metode Adams
        std::vector<double> next_y = adams_moulton_step(time[i], q_charge[i], i_current[i], h, f_history);
        q_charge[i+1] = next_y[0];
        i_current[i+1] = next_y[1];
        
        // Update histori turunan: hapus yang paling lama, tambahkan yang baru
        f_history.erase(f_history.begin());
        f_history.push_back({f1(time[i+1], q_charge[i+1], i_current[i+1]), f2(time[i+1], q_charge[i+1], i_current[i+1])});
    }

    // --- MENULIS HASIL KE FILE ---
    // Membuat file output .csv untuk kemudahan plotting dan analisis
    std::ofstream output_file("adams_moulton_results.csv");
    output_file << "Time,Charge,Current" << std::endl;
    for (int i = 0; i <= n_steps; ++i) {
        output_file << time[i] << "," << q_charge[i] << "," << i_current[i] << std::endl;
    }
    output_file.close();

    std::cout << "Simulasi selesai. Hasil disimpan di adams_moulton_results.csv" << std::endl;

    return 0;
}