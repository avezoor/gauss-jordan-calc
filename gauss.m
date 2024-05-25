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
  m = input('Masukkan jumlah baris (m): ');
  n = input('Masukkan jumlah kolom (n): ');
  fprintf('Masukkan elemen-elemen matriks:\n');
  A = zeros(m, n);
  for i = 1:m
    for j = 1:n
      A(i, j) = input(sprintf('Masukkan elemen A(%d, %d): ', i, j));
    end
  end
  clc
  fprintf('Author: <OBERON AVEZOOR>\n');
  fprintf('Project: Gauss Jordan Elimination Calculator\n');
  fprintf('Created: 2024-05-26>\n');
  fprintf('Augmented Matrix yang diinputkan:\n');
  fprintf('Augmented Matriks :\n');
  disp(A);
  fprintf('\n');
  for i = 1:min(m, n-1)
  if A(i, i) == 0
      [~, max_idx] = max(abs(A(i:m, i)));
      max_idx = max_idx + i - 1;
      if max_idx ~= i
        A([i, max_idx], :) = A([max_idx, i], :);
        fprintf('Tukar B%d dengan B%d:\n', i, max_idx);
          disp(A);
        fprintf('\n');
      end
    end
    if A(i, i) ~= 0 && A(i, i) ~= 1
      multiplier = A(i, i);
        A(i, :) = A(i, :) / A(i, i);
          fprintf('Menjadikan B%d memiliki L1 dengan multiplier %s:\n', i, to_fraction(multiplier));
        disp(A);
      fprintf('\n');
    end
    for j = 1:m
      if j ~= i && A(j, i) ~= 0
        multiplier = A(j, i);
        if i == m && all(A(i, 1:n-1) == 0)
          fprintf('Baris terakhir membentuk elemen nol semua. Proses eliminasi dihentikan.\n');
          fprintf('Solusi variabel:\n');
          K = zeros(1, n-1);
          for i = 1:n-1
            if i <= m
              K(i) = A(i, end);
              fprintf('K%d = %s\n', i, to_fraction(K(i)));
            else
              fprintf('K%d = 0\n', i);
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
            A(j, :) = A(j, :) - A(j, i) * A(i, :);
          disp(A);
          fprintf('\n');
        end
      end
    end
    if i == m && all(A(i, 1:n-1) == 0)
      fprintf('Baris terakhir membentuk elemen nol semua. Proses eliminasi dihentikan.\n');
      fprintf('Solusi variabel:\n');
      K = zeros(1, n-1);
      for i = 1:n-1
        if i <= m
          K(i) = A(i, end);
          fprintf('K%d = %s\n', i, to_fraction(K(i)));
        else
          fprintf('K%d = 0\n', i);
        end
      end
      fprintf('\n');
      return;
    end
  end
  if A(m, n) ~= 0
    fprintf('Sistem tidak konsisten\n');
    disp(A);
    fprintf('Solusi variabel:\n');
    K = zeros(1, n-1);
    for i = 1:n-1
      if i <= m
        K(i) = A(i, end);
        fprintf('K%d = %s\n', i, to_fraction(K(i)));
      else
        fprintf('K%d = 0\n', i);
      end
    end
    free_rows = [];
    for i = 1:m
      if all(A(i, 1:n-1) == 0) && A(i, n) ~= 0
        free_rows = [free_rows, i];
      end
    end
    if ~isempty(free_rows)
      fprintf('Baris free: %s\n', num2str(free_rows));
      fprintf('Proses eliminasi dihentikan.\n');
      return;
    end
    fprintf('\n');
    return;
  end
  if all(A(m, 1:n-1) == 0)
    fprintf('Semua elemen kolom terakhir pada baris terakhir adalah nol. Proses eliminasi dihentikan.\n');
    fprintf('Solusi variabel:\n');
    K = zeros(1, n-1);
    for i = 1:n-1
      if i <= m
        K(i) = A(i, end);
        fprintf('K%d = %s\n', i, to_fraction(K(i)));
      else
        fprintf('K%d = 0\n', i);
      end
    end
    return;
  end
end

function fraction_str = to_fraction(num)
  [numerator, denominator] = rat(num);
  if denominator == 1
    fraction_str = sprintf('%d', numerator);
  else
    fraction_str = sprintf('%d/%d', numerator, denominator);
  end
end

gauss();
