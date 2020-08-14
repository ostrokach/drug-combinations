# The parallel way
parallel -j 8 python chemical_interactions_v2.py --output_db /media/alexey/ssd/gridsearch.db --path_to_data /home/kimlab1/strokach/databases/chemical_interactions/ --input_file version_2.1/predictors_2/predictor_1.tsv --clf_type 1 --n_folds 40 --n_estimators {1} --learning_rate {2}  --max_depth {3} --subsample {4} ::: 25 50 100 200 300 400 500 600 750 1000 1250 1500 1750 2000 2500 3000 ::: 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.1 0.12 0.15 0.2 0.4 0.5 0.6 0.8 1.0 ::: 2 3 4 5 6 ::: 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0


# The loopy way
do for n_estimators in 25 50 100 200 300 400 500 600 750 1000 1250 1500 1750 2000 2500 3000
	do for learning_rate in 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.10 0.12 0.15 0.20 0.40 0.50 0.60 0.80 1.00 1.2
		do for max_depth in 2 3 4 5 6
			do for subsample in 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
				do 
				sleep 1
				done
			done
		done
	done
done


#

