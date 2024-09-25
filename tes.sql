SELECT * FROM hrdataset_lintang
-- Case Study --> Business Queastion
-- 1. Berapa rata-rata salary dari setiap karyawan per department?
-- 2. Berapa banyak karyawan yang sudah bekerja lebih dari 8 tahun?
-- 3. Ada korelasi antara performance score karyawan dengan level kepuasan?
-- 4. Apakah ada hubungan antara sering absen dengan performance score?
-- No 1
SELECT department, AVG(salary) as avg_salary
FROM hrdataset_lintang
GROUP BY department
ORDER BY avg_salary desc;
-- No 2
SELECT COUNT (*) AS number_of_employees
FROM hrdataset_lintang
WHERE DATE_PART('year', AGE(dateofhire::date)) > 8
-- No 3
SELECT performancescore, AVG(empsatisfaction) as Avg_satisfaction
FROM hrdataset_lintang
GROUP BY performancescore
ORDER BY avg_satisfaction;
-- Kesimpulan: Berkorelasi positif(yg satu naik, satunya juga naik), korelasi negatif kebalikannya (satu naik, yg satunya turun). Kalo employee nya bagus performancenya (exceeds), nilainya juga tinggi
-- PIP: employee masih perlu improvement
-- No 4
SELECT absences, AVG(perfscoreid) as avg_performance
FROM hrdataset_lintang
GROUP BY absences
ORDER BY absences;
-- Kesimpulan: tidak ada korelasi yang significant