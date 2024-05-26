## Copyright (C) 2024 Oberon Avezoor
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.
## Author: oberon <oberon@AVEZOOR>
## Created: 2024-05-26
clc
#-------------------------------------------------------------------------------
  % Meminta bentuk Pecahan
#------------------------------------------------------------------------------
format rat
tic;
function gauss()
#-------------------------------------------------------------------------------
  % Meminta input ukuran matriks m x n dari pengguna
#-------------------------------------------------------------------------------
  m = input('Masukkan jumlah baris (m): ');
  n = input('Masukkan jumlah kolom (n): ');
#-------------------------------------------------------------------------------
  % Meminta input elemen-elemen matriks dari pengguna dan definisi langsung
#-------------------------------------------------------------------------------
  fprintf('Masukkan elemen-elemen matriks:\n');
  A = zeros(m, n);
  for i = 1:m
    for j = 1:n
      A(i, j) = input(sprintf('Masukkan elemen A(%d, %d): ', i, j));
    end
  end
#-------------------------------------------------------------------------------
  % Menampilkan matriks awal atau Augmented Matriks
#-------------------------------------------------------------------------------
  clc
  fprintf('---------------------------------------------------\n');
  fprintf('|  Author: <oberon@AVEZOOR>                       |\n');
  fprintf('|  Project: Gauss Jordan Elimination Calculator   |\n');
  fprintf('|  Created: 2024-05-26>                           |\n');
  fprintf('---------------------------------------------------\n');
  fprintf('Augmented Matrix yang diinputkan:\n');
  disp(A);
  fprintf('\n');
#-------------------------------------------------------------------------------
  % Iterasi untuk setiap kolom hingga kolom terakhir yang relevan (n-1)
  % Partial pivoting jika elemen diagonal adalah nol
  % Mencari nilai maksimum dalam kolom untuk pivot
  % Mengubah indeks relatif ke indeks absolut
  % Menukar baris jika pada kolom lending one bernilai nol
#-------------------------------------------------------------------------------
  for i = 1:min(m, n-1)
    if A(i, i) == 0
      [~, max_idx] = max(abs(A(i:m, i)));
      max_idx = max_idx + i - 1;
      if max_idx ~= i
        A([i, max_idx], :) = A([max_idx, i], :);
        fprintf('Tukar Baris B%d dengan B%d:\n', i, max_idx);
        disp(A);
        fprintf('\n');
      end
    end
#-------------------------------------------------------------------------------
  % Membuat elemen diagonal menjadi 1
  % Menyimpan nilai diagonal sebelum dibagi
  % Membagi seluruh baris untuk membuat diagonal menjadi 1
#-------------------------------------------------------------------------------
    if A(i, i) ~= 0 && A(i, i) ~= 1
      multiplier = A(i, i);
      A(i, :) = A(i, :) / A(i, i);
      fprintf('Menjadikan L1 pada B%d | Multiplier %s\n', i, to_fraction(multiplier));
      disp(A);
      fprintf('\n');
    end
#-------------------------------------------------------------------------------
  % Membuat elemen lainnya di kolom tersebut menjadi 0
  % Menyimpan nilai elemen yang akan di-nolkan
  % Hentikan looping OBE jika baris terakhir membentuk elemen nol semua
  % Menampilkan solusi variabel K1, K2, K3, ...,Kn
  % Larik untuk menyimpan nilai variabel K
  % Menampilkan solusi dalam bentuk pecahan
  % Jika tidak ada solusi (baris free), set ke free
#-------------------------------------------------------------------------------
    for j = 1:m
      if j ~= i && A(j, i) ~= 0
        multiplier = A(j, i);
        if i == m && all(A(i, 1:n-1) == 0)
          fprintf('Proses eliminasi dihentikan. Terima kasih\n');
          fprintf('---------------------------------------------------\n');
          fprintf('Solusi variabel:\n');
          K = zeros(1, n-1);
          for i = 1:n-1
            if i <= m
              K(i) = A(i, end);
              fprintf('K%d = %s\n', i, to_fraction(K(i)));
            else
              fprintf('K%d = free\n', i);
            end
          end
          fprintf('\n');
          return;
        end
        if multiplier ~= 0
          multiplier_str = to_fraction(abs(multiplier));
          if multiplier < 0
            fprintf('OBE pada B%d%d(%s)\n', j, i, multiplier_str);
          else
            fprintf('OBE pada B%d%d(-%s)\n', j, i, multiplier_str);
          end
          A(j, :) = A(j, :) - A(j, i) * A(i, :);
          disp(A);
          fprintf('\n');
        end
      end
    end
#-------------------------------------------------------------------------------
  % Hentikan looping OBE jika baris terakhir membentuk elemen nol semua
  % Menampilkan solusi variabel K1, K2, K3, ...,Kn
  % Larik untuk menyimpan nilai variabel K
  % Menampilkan solusi dalam bentuk pecahan
  % Jika tidak ada solusi (baris free), set ke free
#-------------------------------------------------------------------------------
    if i == m && all(A(i, 1:n-1) == 0)
      fprintf('Proses eliminasi dihentikan. Terima kasih\n');
      fprintf('---------------------------------------------------\n');
      fprintf('Solusi variabel:\n');
      K = zeros(1, n-1);
      for i = 1:n-1
        if i <= m
          K(i) = A(i, end);
          fprintf('K%d = %s\n', i, to_fraction(K(i)));
        else
          fprintf('K%d = free\n', i);
        end
      end
      fprintf('\n');
      return;
    end
  end
#-------------------------------------------------------------------------------
  % Hentikan proses jika kolom terakhir pada baris akhir adalah angka bukan nol
  % Menampilkan solusi variabel K1, K2, K3, ...,Kn
  % Larik untuk menyimpan nilai variabel K
  % Menampilkan solusi dalam bentuk pecahan
#-------------------------------------------------------------------------------
  if A(m, n) ~= 0
    fprintf('Sistem tidak konsisten\n');
    disp(A);
    fprintf('---------------------------------------------------\n');
    fprintf('Solusi variabel:\n');
    K = zeros(1, n-1);
    for i = 1:n-1
      if i <= m
        K(i) = A(i, end);
        fprintf('K%d = %s\n', i, to_fraction(K(i)));
      else
        fprintf('K%d = free\n', i);
      end
    end
#-------------------------------------------------------------------------------
    % Cari baris yang free (semua elemen di kolom utama adalah nol, namun hasilnya tidak nol)
#-------------------------------------------------------------------------------
    free_rows = [];
    for i = 1:m
      if all(A(i, 1:n-1) == 0) && A(i, n) ~= 0
        free_rows = [free_rows, i];
      end
    end
#-------------------------------------------------------------------------------
    % Hentikan proses jika sistem menjadi inkonsisten
#-------------------------------------------------------------------------------
    if ~isempty(free_rows)
      fprintf('---------------------------------------------------\n');
      fprintf('Baris free: %s\n', num2str(free_rows));
      fprintf('Proses eliminasi dihentikan. Terima kasih\n');
      return;
    end
    fprintf('\n');
    return;
  end
#-------------------------------------------------------------------------------
  % Hentikan proses jika semua elemen kolom terakhir pada baris akhir adalah nol
#-------------------------------------------------------------------------------
  if all(A(m, 1:n-1) == 0)
    fprintf('Proses eliminasi dihentikan. Terima kasih\n');
    fprintf('---------------------------------------------------\n');
    fprintf('Solusi variabel:\n');
    K = zeros(1, n-1);
    for i = 1:n-1
      if i <= m
        K(i) = A(i, end);
        fprintf('K%d = %s\n', i, to_fraction(K(i)));
      else
        fprintf('K%d = free\n', i);
      end
    end
    return;
  end
end
#-------------------------------------------------------------------------------
  % Fungsi untuk mengonversi pecahan menjadi n/m
  % Mencari penyebut terkecil untuk mendapatkan pembilang dan penyebut yang sederhana
  % Hanya menampilkan pembilang jika penyebut adalah 1
#-------------------------------------------------------------------------------
function fraction_str = to_fraction(num)
  [numerator, denominator] = rat(num);
  if denominator == 1
    fraction_str = sprintf('%d', numerator);
  else
    fraction_str = sprintf('%d/%d', numerator, denominator);
  end
end
#-------------------------------------------------------------------------------
% Panggil fungsi gauss
#-------------------------------------------------------------------------------
gauss();
excuted_time = toc;
fprintf('---------------------------------------------------\n');
fprintf('Waktu eksekusi: %.6f detik\n', excuted_time);
