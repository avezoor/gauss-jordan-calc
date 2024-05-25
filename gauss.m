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
  fprintf('Author: <oberon@AVEZOOR>\n');
  fprintf('Project: Gauss Jordan Elimination Calculator\n');
  fprintf('Created: 2024-05-26>\n');
  fprintf('Augmented Matrix yang diinputkan:\n');
  augmatrix(A);
  fprintf('\n');

  % Iterasi untuk setiap kolom hingga kolom terakhir yang relevan (n-1)
  for i = 1:min(m, n-1)
    % Partial pivoting jika elemen diagonal adalah nol
    if A(i, i) == 0
      [~, maxindex] = max(abs(A(i:m, i))); % Mencari nilai maksimum dalam kolom untuk pivot
      maxindex = maxindex + i - 1; % Mengubah indeks relatif ke indeks absolut
      if maxindex ~= i
        A([i, maxindex], :) = A([maxindex, i], :); % Menukar baris
        fprintf('Tukar B%d dengan B%d:\n', i, maxindex);
        augmatrix(A);
        fprintf('\n');
      end
    end

    % Membuat elemen diagonal menjadi 1 jika belum
    if A(i, i) ~= 0 && A(i, i) ~= 1
      multiplier = A(i, i); % Menyimpan nilai diagonal sebelum dibagi
      A(i, :) = A(i, :) / A(i, i); % Membagi seluruh baris untuk membuat diagonal menjadi 1
      fprintf('Menjadikan B%d memiliki L1 dengan multiplier %s:\n', i, tofrac(multiplier));
      augmatrix(A);
      fprintf('\n');
    end

    % Membuat elemen lainnya di kolom tersebut menjadi 0
    for j = 1:m
      if j ~= i && A(j, i) ~= 0
        multiplier = A(j, i); % Menyimpan nilai elemen yang akan di-nolkan
        if multiplier ~= 0
          fprintf('OBE pada B%d %s (%s)B%d\n', j, plusmin(multiplier), tofrac(abs(multiplier)), i);
          A(j, :) = A(j, :) - A(j, i) * A(i, :); % Mengurangi baris untuk membuat elemen menjadi 0
          augmatrix(A);
          fprintf('\n');
        end
      end
    end
  end

  % Cari baris yang free (semua elemen di kolom utama adalah nol, namun hasilnya tidak nol)
  isfree = [];
  for i = 1:m
    if all(A(i, 1:n-1) == 0) && A(i, n) ~= 0
      isfree = [isfree, i];
    end
  end

  % Menampilkan solusi variabel K1, K2, K3, ...
  fprintf('Solusi variabel:\n');
  for i = 1:n-1
    if i <= m
      fprintf('K%d = %s\n', i, tofrac(A(i, end))); % Menampilkan solusi dalam bentuk pecahan
    else
      fprintf('K%d = 0\n', i); % Jika tidak ada solusi (baris free), set ke 0
    end
  end
  fprintf('\n');

  % Hentikan proses jika sistem menjadi inkonsisten
  if ~isempty(isfree)
    fprintf('Sistem menjadi inkonsisten: Baris free: %s\n', num2str(isfree));
    fprintf('Proses eliminasi dihentikan.\n');
    return;
  end

  % Menampilkan matriks akhir dalam bentuk Reduced Row Echelon Form (RREF)
  fprintf('Matriks akhir dalam bentuk Reduced Row Echelon Form (RREF):\n');
  augmatrix(A);
end

% Fungsi untuk menampilkan matriks
function augmatrix(A)
  disp(A);
end

% Fungsi untuk mengonversi pecahan menjadi n/m
function fraction_str = tofrac(num)
  % Mencari penyebut terkecil untuk mendapatkan pembilang dan penyebut yang sederhana
  [numerator, denominator] = rat(num);
  if denominator == 1
    fraction_str = sprintf('%d', numerator); % Hanya menampilkan pembilang jika penyebut adalah 1
  else
    fraction_str = sprintf('%d/%d', numerator, denominator);
  end
end

% Fungsi untuk menampilkan tanda +/- sesuai dengan nilai multiplier
function str = plusmin(multiplier)
  if multiplier < 0
    str = '-';
  else
    str = '+';
  end
end

% Panggil fungsi gauss
gauss();
