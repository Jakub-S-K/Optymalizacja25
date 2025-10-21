#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wej�ciowe:
	// ff - wska�nik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i g�rne ograniczenie
	// epslion - zak��dana dok�adno�� rozwi�zania
	// Nmax - maksymalna liczba wywo�a� funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy warto�� funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi�zanie z zadan� dok�adno�ci�
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo�a� funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        double* p = new double[2] { 0, 0 };
 
        int i = 0;
        solution Xopt0 = x0;
        solution Xopt1(x0 + d);
 
        Xopt0.fit_fun(ff, ud1, ud2);
        Xopt1.fit_fun(ff, ud1, ud2);
 
        if (Xopt1.y == Xopt0.y)
        {
            p[0] = m2d(Xopt0.x);
            p[1] = m2d(Xopt1.x);
            return p;
        }
 
        if (Xopt1.y > Xopt0.y)
        {
            d = -d;
            Xopt1.x = Xopt0.x + d;
 
            Xopt1.fit_fun(ff, ud1, ud2);
            if (Xopt1.y >= Xopt0.y)
            {
                p[0] = m2d(Xopt1.x);
                p[1] = m2d(Xopt0.x - d);
                return p;
            }
        }
 
        solution Xopt_prev;
        while (true)
        {
            if (solution::f_calls > Nmax)
            {
                Xopt0.flag = 0;                                 // flaga = 0 ozancza przekroczenie maksymalne liczby wywołań funkcji celu
                break;
            }
 
            Xopt_prev = Xopt1;
 
            i++;
            Xopt1.x = Xopt0.x + pow(alpha, i) * d;
            Xopt1.fit_fun(ff, ud1, ud2);
 
            if (Xopt1.y > Xopt_prev.y) //na ODWROT???? NIE, DZIAŁA JEDNAK!!!!!! OSZALEJE HEHAHAHAHEHAFADOS;IFHA;SDKJ AKSDLASDF
                break;
        }
 
        if (d > 0)
        {
            p[0] = m2d(Xopt_prev.x);
            p[1] = m2d(Xopt1.x);
        }
        else
        {
            p[0] = m2d(Xopt1.x);
            p[1] = m2d(Xopt_prev.x);
        }
 
        return p;
    }
    catch (string ex_info)
    {
        throw ("double* expansion(...):\n" + ex_info);
    }
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* interval;
		interval = expansion(ff, a, 5, 1.2, Nmax);
		solution Xopt;
		solution::clear_calls();
		double aa[3], bb[3], cc[3], dd[3];
		int i = 1;
		aa[i] = interval[0];
		bb[i] = interval[1];
		cc[i] = (interval[1] + interval[0]) / 2;

		do {
			Xopt.x = aa[i];
			Xopt.fit_fun(ff, ud1, ud2);
			double l_f_a = m2d(Xopt.y);

			Xopt.x = bb[i];
			Xopt.fit_fun(ff, ud1, ud2);
			double l_f_b = m2d(Xopt.y);
			
			Xopt.x = cc[i];
			Xopt.fit_fun(ff, ud1, ud2);
			double l_f_c = m2d(Xopt.y);

			double l = l_f_a * (pow(bb[i], 2) - pow(cc[i], 2)) + l_f_b * (pow(cc[i], 2) - pow(aa[i], 2)) + l_f_c * (pow(aa[i], 2) - pow(bb[i], 2));

			double m = l_f_a * (bb[i] - cc[i]) + l_f_b * (cc[i] - aa[i]) + l_f_c * (aa[i] - bb[i]);
			if (m <= 0) {
				throw std::string("M <= 0");
			}
			
			dd[i] = 0.5 * l / m;

			Xopt.x = dd[i];
			Xopt.fit_fun(ff, ud1, ud2);
			double l_f_d = m2d(Xopt.y);

			if ((aa[i] < dd[i]) && (dd[i] < cc[i])) {
				if (l_f_d < l_f_c) {
					aa[i + 1] = aa[i];
					cc[i + 1] = dd[i];
					bb[i + 1] = cc[i];
				} else {
					aa[i + 1] = dd[i];
					cc[i + 1] = cc[i];
					bb[i + 1] = bb[i];
				}
			} else {
				if ((cc[i] < dd[i]) && (dd[i] < bb[i])) {
					if (l_f_d < l_f_c) {
						aa[i + 1] = cc[i];
						cc[i + 1] = dd[i];
						bb[i + 1] = bb[i];
					} else {
						aa[i + 1] = aa[i];
						cc[i + 1] = cc[i];
						bb[i + 1] = dd[i];
					}
				} else {
					throw std::string("Error 4321");
				}
			}
			
			aa[i - 1] = aa[i];
			bb[i - 1] = bb[i];
			cc[i - 1] = cc[i];
			dd[i - 1] = dd[i];
						
			aa[i] = aa[i + 1];
			bb[i] = bb[i + 1];
			cc[i] = cc[i + 1];
			dd[i] = dd[i + 1];

			if (solution::f_calls > Nmax) {
				throw std::string("Too many interations");
			}

		} while (bb[i] - aa[i] > epsilon && abs(dd[i] - dd[i - 1]) > gamma);

		std::cout << "Wynik: " << dd[i - 1];
		// std::cout << "Test - x: " << Xopt.x << " y: " << Xopt.y << "\n";
		
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
