# Adams-Moulton RLC Simulation

**Nama:** Shafwan Hasyim

**NPM:** 2306209113

---

## Deskripsi Proyek

Proyek ini merupakan implementasi dan analisis metode numerik untuk menyelesaikan persamaan diferensial biasa (ODE) yang merepresentasikan respons transien dari sebuah rangkaian RLC seri. 

Fokus utama dari studi ini adalah untuk menganalisis secara mendalam performa, khususnya karakteristik stabilitas dan akurasi, dari **metode predictor-corrector Adams-Moulton orde empat**. Sebagai studi kasus, digunakan model rangkaian RLC ideal tanpa redaman (R=0) dengan eksitasi sinusoidal, sebuah sistem osilasi non-disipatif yang dapat secara efektif menguji batas kemampuan sebuah algoritma numerik. 

Hasil dari metode Adams-Moulton divalidasi terhadap solusi analitik yang eksak dan luaran dari metode Runge-Kutta Orde Empat (RK4) untuk mengukur performanya secara kuantitatif. 

## Isi Repositori

Repositori ini berisi program utama:

* `main.cpp`: Implementasi utama dalam C++ yang menggunakan metode **Adams-Moulton orde empat**. 

## Cara Menggunakan

Kedua program dapat dikompilasi dan dijalankan menggunakan compiler C++ standar (G++).

**Program main Adams-Moulton (C++)**

```bash
# Kompilasi
g++ main.cpp -o main

# Jalankan
./main
```

## Output Program

Setelah program berhasil dijalankan, sebuah file baru bernama adams_moulton_results.csv akan dibuat di dalam direktori yang sama.

File ini berformat CSV (Comma-Separated Values) yang dapat dibuka dengan mudah menggunakan aplikasi spreadsheet seperti Microsoft Excel atau Google Sheets. File ini berisi tiga kolom data:

1. Time: Waktu simulasi (dalam detik).
2. Charge: Nilai muatan $q(t)$ pada kapasitor (dalam Coulomb).
3. Current: Nilai arus $i(t)$ yang mengalir dalam rangkaian (dalam Ampere).
