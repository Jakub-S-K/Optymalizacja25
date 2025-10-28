#include"user_funs.h"
#include <math.h>

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsprzdne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocztkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment siy dziaajcy na wahado oraz czas dziaania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwizujemy rwnanie rniczkowe
	int n = get_len(Y[0]);									// dugo rozwizania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahada
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto funkcji celu (ud1 to zaoone maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pamici rozwizanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z pooenia to prdko
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z prdkoci to przyspieszenie
	return dY;
}

matrix ff_test(matrix x_m, matrix ud1, matrix ud2) {
	double x = m2d(x_m);
	double exponent = - ((0.1 * x - 2*M_PI)*(0.1 * x - 2*M_PI));
	x = 0 - cos(0.1 * x) * pow(M_E, exponent) + 0.002 * (0.1 * x) * (0.1 * x);
	return matrix(x);
}

matrix df1(double x, matrix Y, matrix ud1, matrix ud2)
{
    // Y = [Va, Vb, Tb]
    matrix dY(3, 1);
    double a = 0.98, b = 0.63, g = 9.81;
    double Pa = 2.0, Pb = 1.0;
    double Ta = 95.0;                 // temperatura w zbiorniku A
    double Tin = 20.0;                // temperatura dopływu z zewnątrz
    double Fin = 0.01;                // dopływ z zewnątrz [m^3/s]
    double Db = 0.00365665;           // [m^2]

    double Da = m2d(ud1);             // [m^2] pole otworu między A a B

    double Fa_out = a*b*Da*sqrt((2*g*Y(0))/Pa);       // wypływ z A do B
    double Fb_out = a*b*Db*sqrt((2*g*Y(1))/Pb);       // wypływ z B na zewnątrz

    // Równania różniczkowe:
    dY(0) = -Fa_out;
    dY(1) = Fa_out + Fin - Fb_out;
    dY(2) = (Fin/Y(1))*(Tin - Y(2)) + (Fa_out/Y(1))*(Ta - Y(2));

    return dY;
}


matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
    // x - wektor argumentów (x(0) = D_A w cm^2)
    // ud1, ud2 - nieużywane (opcjonalne)
    
    matrix y(1,1);
    matrix Y0 = matrix(3, 1);
    Y0(0) = 5.0;    // Va0 [m^3]
    Y0(1) = 1.0;    // Vb0 [m^3]
    Y0(2) = 20.0;   // Tb0 [°C]
    
    // Rozwiązanie ODE
    matrix* Y = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, x, NULL);

    // Poprawka: długość (liczba kroków) pobierana z wektora czasu Y[0]
    int n = get_len(Y[0]); // liczba wierszy (kroków czasowych)
    if (n <= 0) {
        // jawny komunikat w razie błędnego/niepełnego rozwiązania ODE
        Y[0].~matrix();
        Y[1].~matrix();
        throw string("ff1R: nieprawidlowe rozwiazanie ODE (n <= 0)");
    }
    
    // Szukamy maksymalnej temperatury w zbiorniku B -> trzecia zmienna stanu to kolumna 2
    double Tmax = Y[1](0, 2);
    for (int i = 1; i < n; ++i)
        if (Y[1](i, 2) > Tmax)
            Tmax = Y[1](i, 2);
    
    // Funkcja celu: różnica względem 50°C (ud1 zawiera cel temperatury)
    y(0,0) = fabs(Tmax - m2d(ud1));
    
    // Zwalniamy pamięć macierzy zgodnie z innymi funkcjami w projekcie
    Y[0].~matrix();
    Y[1].~matrix();
    
    return y;
}


matrix ff2T(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));

    return matrix(x1*x1 + x2*x2 - cos(2.5*M_PI*x1) - cos(2.5*M_PI*x2) + 2);
}

// matrix ff2T(matrix x, matrix ud1, matrix ud2) {
//     double x1 = m2d(x(0));
//     double x2 = m2d(x(1));

//     return matrix(x1*x1 + x2*x2 - cos(2.5*M_PI*x1) - cos(2.5*M_PI*x2) + 2);
// }