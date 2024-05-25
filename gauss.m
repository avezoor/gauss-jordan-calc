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

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} g (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: oberon <oberon@AVEZOOR>
## Created: 2024-05-26

function gauss()
  % Meminta input ukuran matriks m x n dari pengguna
  m = input('Masukkan jumlah baris (m): ');
  n = input('Masukkan jumlah kolom (n): ');

  % Meminta input elemen-elemen matriks dari pengguna dan definisi langsung
  fprintf('Masukkan elemen-elemen matriks:\n');
  A = zeros(m, n);
  for i = 1:m
    for j = 1:n
      A(i, j) = input(sprintf('Masukkan elemen A(%d, %d): ', i, j));
    end
  end

  % Menampilkan matriks awal
  clc
  fprintf('Author: <OBERON AVEZOOR>\n');
  fprintf('Project: Gauss Jordan Elimination Calculator\n');
  fprintf('Created: 2024-05-26>\n');
  fprintf('Augmented Matrix yang diinputkan:\n');
  fprintf('Augmented Matriks :\n');
  disp(A);
  fprintf('\n');

  % Iterasi untuk setiap kolom hingga kolom terakhir yang relevan (n-1)
  for i = 1:min(m, n-1)
    % Partial pivoting jika elemen diagonal adalah nol
    if A(i, i) == 0
      [~, max_idx] = max(abs(A(i:m, i))); % Mencari nilai maksimum dalam kolom untuk pivot
      max_idx = max_idx + i - 1; % Mengubah indeks relatif ke indeks absolut
      if max_idx ~= i
        A([i, max_idx], :) = A([max_idx, i], :); % Menukar baris
        fprintf('Tukar B%d dengan B%d:\n', i, max_idx);
        disp(A);
        fprintf('\n');
      end
    end

    % Membuat elemen diagonal menjadi 1 jika belum
    if A(i, i) ~= 0 && A(i, i) ~= 1
      multiplier = A(i, i); % Menyimpan nilai diagonal sebelum dibagi
      A(i, :) = A(i, :) / A(i, i); % Membagi seluruh baris untuk membuat diagonal menjadi 1
      fprintf('Menjadikan B%d memiliki L1 dengan multiplier %s:\n', i, to_fraction(multiplier));
      disp(A);
      fprintf('\n');
    end

    % Membuat elemen lainnya di kolom tersebut menjadi 0
    for j = 1:m
      if j ~= i && A(j, i) ~= 0
        multiplier = A(j, i); % Menyimpan nilai elemen yang akan di-nolkan
            % Hentikan looping OBE jika baris terakhir membentuk elemen nol semua
        if i == m && all(A(i, 1:n-1) == 0)
          fprintf('Baris terakhir membentuk elemen nol semua. Proses eliminasi dihentikan.\n');
          % Menampilkan solusi variabel K1, K2, K3, ...
          fprintf('Solusi variabel:\n');
          K = zeros(1, n-1); % Larik untuk menyimpan nilai variabel K
          for i = 1:n-1
            if i <= m
              K(i) = A(i, end); % Simpan nilai variabel K
              fprintf('K%d = %s\n', i, to_fraction(K(i))); % Menampilkan solusi dalam bentuk pecahan
            else
              fprintf('K%d = 0\n', i); % Jika tidak ada solusi (baris free), set ke 0
            end
          end
          fprintf('\n');
          return;
        end
        if multiplier ~= 0
          multiplier_str = to_fraction(abs(multiplier));
          if multiplier < 0
            fprintf('OBE pada B%d + (%s)B%d:\n', j, multiplier_str, i);
          else
            fprintf('OBE pada B%d - (%s)B%d:\n', j, multiplier_str, i);
          end
          A(j, :) = A(j, :) - A(j, i) * A(i, :); % Mengurangi baris untuk membuat elemen menjadi 0
          disp(A);
          fprintf('\n');
        end
      end
    end

    % Hentikan looping OBE jika baris terakhir membentuk elemen nol semua
    if i == m && all(A(i, 1:n-1) == 0)
      fprintf('Baris terakhir membentuk elemen nol semua. Proses eliminasi dihentikan.\n');
      % Menampilkan solusi variabel K1, K2, K3, ...
      fprintf('Solusi variabel:\n');
      K = zeros(1, n-1); % Larik untuk menyimpan nilai variabel K
      for i = 1:n-1
        if i <= m
          K(i) = A(i, end); % Simpan nilai variabel K
          fprintf('K%d = %s\n', i, to_fraction(K(i))); % Menampilkan solusi dalam bentuk pecahan
        else
          fprintf('K%d = 0\n', i); % Jika tidak ada solusi (baris free), set ke 0
        end
      end
      fprintf('\n');
      return;
    end
  end

  % Hentikan proses jika kolom terakhir pada baris akhir adalah angka bukan nol
  if A(m, n) ~= 0
    fprintf('Sistem tidak konsisten\n');
    disp(A);
    % Menampilkan solusi variabel K1, K2, K3, ...
    fprintf('Solusi variabel:\n');
    K = zeros(1, n-1); % Larik untuk menyimpan nilai variabel K
    for i = 1:n-1
      if i <= m
        K(i) = A(i, end); % Simpan nilai variabel K
        fprintf('K%d = %s\n', i, to_fraction(K(i))); % Menampilkan solusi dalam bentuk pecahan
      else
        fprintf('K%d = 0\n', i); % Jika tidak ada solusi (baris free), set ke 0
      end
    end
    % Cari baris yang free (semua elemen di kolom utama adalah nol, namun hasilnya tidak nol)
    free_rows = [];
    for i = 1:m
      if all(A(i, 1:n-1) == 0) && A(i, n) ~= 0
        free_rows = [free_rows, i];
      end
    end

    % Hentikan proses jika sistem menjadi inkonsisten
    if ~isempty(free_rows)
      fprintf('Baris free: %s\n', num2str(free_rows));
      fprintf('Proses eliminasi dihentikan.\n');
      return;
    end
    fprintf('\n');
    return;
  end

  % Hentikan proses jika semua elemen kolom terakhir pada baris akhir adalah nol
  if all(A(m, 1:n-1) == 0)
    fprintf('Semua elemen kolom terakhir pada baris terakhir adalah nol. Proses eliminasi dihentikan.\n');
        % Menampilkan solusi variabel K1, K2, K3, ...
    fprintf('Solusi variabel:\n');
    K = zeros(1, n-1); % Larik untuk menyimpan nilai variabel K
    for i = 1:n-1
      if i <= m
        K(i) = A(i, end); % Simpan nilai variabel K
        fprintf('K%d = %s\n', i, to_fraction(K(i))); % Menampilkan solusi dalam bentuk pecahan
      else
        fprintf('K%d = 0\n', i); % Jika tidak ada solusi (baris free), set ke 0
      end
    end
    return;
  end

end

% Fungsi untuk mengonversi pecahan menjadi n/m
function fraction_str = to_fraction(num)
  % Mencari penyebut terkecil untuk mendapatkan pembilang dan penyebut yang sederhana
  [numerator, denominator] = rat(num);
  if denominator == 1
    fraction_str = sprintf('%d', numerator); % Hanya menampilkan pembilang jika penyebut adalah 1
  else
    fraction_str = sprintf('%d/%d', numerator, denominator);
  end
end

% Panggil fungsi gauss
gauss();
