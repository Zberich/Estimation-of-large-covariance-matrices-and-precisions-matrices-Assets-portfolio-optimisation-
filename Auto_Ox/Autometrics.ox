#include <oxstd.oxh>
#import <packages/PcGive/pcgive>

// Function running Autometrics with "classic" standard error
run_Autometrics_U(const obj, const print_output, const name_Y, const a_names_X, const a_names_U, const Auto, const TS_Auto, const IIS, const TS_IIS,const d_var_name){
	//Add const"ant to the database
	obj.Deterministic(-1);
	obj.DeSelect();
	obj.Select("Y", {name_Y, 0, 0});
	for (decl i=0;i<sizerc(a_names_U);++i){
	 	obj.Select("U", {a_names_U[i], 0, 0});
	}

	for (decl i=0;i<sizerc(a_names_X);++i) {
	   	//X variables are variables not in 
		if(find(a_names_U,a_names_X[i])<0)
			obj.Select("X", {a_names_X[i], 0, 0});
	}
	if (Auto){
		obj.Autometrics(TS_Auto, "none", 1);
		obj.AutometricsSet("test_normality", 1);
		obj.AutometricsSet("test_hetero", 1);
		obj.AutometricsSet("test_heterox", 0);
		obj.AutometricsSet("test_chow", 0);
		obj.AutometricsSet("test_reset", 1);
		obj.AutometricsSet("test_ar", 0);
		obj.AutometricsSet("test_portmanteau", 0);
		obj.AutometricsSet("test_arch", 0);
		obj.AutometricsSet("test_default", 0);
	}
	if (IIS)
		obj.Autometrics(TS_IIS, "IIS", 1);
	obj.SetSelSample(-1, 1, -1, 1);
	
	//Estimate  by OLS
	obj.SetMethod("OLS");
	obj.SetPrint(print_output);
	
	obj.Estimate();

	decl R2 = 1- sumsqrc(obj.GetResiduals())/(obj.GetcT()*varc(obj.GetY()));
	
	//Get index of all the selected variables
	decl all_kept_vars_index = 	obj.GetVarIndex(obj.GetParNames()) ;

	//Get index of the U matrix variables
	decl a_index_U = obj.GetVarIndex(a_names_U);

	// kept values are those in the final model and that are not fixed 
	decl kept_vars= exclusion(all_kept_vars_index,a_index_U);

	decl d_var_index = find(obj.GetParNames(),d_var_name) ;

	//Get the variables index, the firsts are Y, we remove -1 to resend indexes in 1 basis form (for R)
	decl index_selvars= kept_vars  	;

	if(d_var_name!="Constant")
		index_selvars= kept_vars-1	 ;


	decl results = obj.GetPi() | index_selvars;
	// In the first row is stored the R2 of the model, the estimated coefficient of d and its standard error,
	// but the index_selvars does not include the idx of D as it is fixed (it only include index of vars selected and not fixed)
		
	return results';
	
}

// Main function reading the elements from R, and applying Autometrics on the input file
main()
{
	decl time_0=timer();
	decl input_file_path="";
	decl output_file_path="";
	decl n_xvars=0;
	decl Auto=1;
	decl TS_Auto=0;
	decl IIS=0;
	decl TS_IIS=0;
	decl idx_alpha=1;
	decl d_var_name = "Constant";
	decl print_output=1;
	
	decl args = arglist(), cargs = sizeof(args);

	if (cargs>1)
	{
		println("args:", args);
		decl idx=1;
        sscan(args[idx++], "%s", &input_file_path);
        sscan(args[idx++], "%d", &n_xvars);
        sscan(args[idx++], "%f", &TS_Auto);
        sscan(args[idx++], "%s", &output_file_path); 
	}
	
	decl model = new PcGive();
	model.Load(input_file_path);
	decl a_names=model.GetAllNames();
	decl name_Y=a_names[0];
	decl name_X_LEVEL=a_names[1:1+n_xvars-1];
	decl name_U_LEVEL={"Constant"};

	decl index_sel=run_Autometrics_U(model,print_output,name_Y,name_X_LEVEL,name_U_LEVEL,Auto,TS_Auto,IIS,TS_IIS,d_var_name);
		savemat(output_file_path,index_sel,{"alpha_hat","vars"});

	print("\nTime elapsed: ", timespan(time_0),"\n"); 
	delete model;	
}

