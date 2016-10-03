#ifndef __RANGEN_H_
#define __RANGEN_H_


class RanGen {
	public:
		//constructor
		RanGen();
		~RanGen();
		//methods
		bool ranbool();
		double randouble();	
		int ranval(int min, int max);	
		double rangauss(double);
		int ranpoisson(double);
		double ranpoissonpdf(unsigned int, double);
		double rancauchypdf(unsigned int, double);
};

#endif
