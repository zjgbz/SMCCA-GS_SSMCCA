module load r/local_3.6.0

for noise_coef in 1 100 100000
do
	for z_coef in 1 100
	do
		for beta_coef in 2 100
		do
			Rscript simple_test_supmcca.R ${noise_coef} ${z_coef} ${beta_coef} > result_${noise_coef}_${z_coef}_${beta_coef}.out 2>&1
		done
	done
done