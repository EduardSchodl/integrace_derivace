
using System;

/**
 * Hlavni trida programu pro numericky vypocet urciteho integralu
 */
public class Exercise04
{
    public static void Main(String[] args)
    {
        /*
        LinearFunction lf = new LinearFunction(5, 1);
        Integrator integrator = new Integrator();
        integrator.SetDelta(0.000001);
        double integralLin = integrator.Integrate(lf, 0, 10);
        Console.WriteLine(integralLin);
        QuadraticPolynomial qp = new QuadraticPolynomial(2, 0, 0);
        double integralQuad = integrator.Integrate(qp, 0, 10);
        Console.WriteLine(integralQuad);

        Console.WriteLine(lf.Differentiate(5));
        Console.WriteLine(qp.Differentiate(5));

        GeneralPolynomial gp = new GeneralPolynomial([7,3,-2,8]);
        Console.WriteLine(integrator.Integrate(gp, 0, 10));
        */

        IFunction myPolynomial = new GeneralPolynomial(new double[] { 7, -5, 3, -15 });
        IFunction firstDerivative = new Derivative(myPolynomial);
        IFunction secondDerivative = new Derivative(firstDerivative);
        IFunction thirdDerivative = new Derivative(secondDerivative);

        Console.WriteLine("třetí:");
        for (double x = 0; x < 10.1; x += 0.2)
        {
            Console.WriteLine(thirdDerivative.ValueAt(x));
        }
        Console.WriteLine("-------------------------");

        // jeste overeni treba pro funkci sinus
        IFunction sin = new Sine();
        IFunction cos = new Derivative(sin);
        Console.WriteLine("Cosinus");
        for (double x = 0; x < Math.PI * 2; x += 0.3)
        {
            Console.WriteLine(cos.ValueAt(x));
        }
    }
}

/**
 * Obecna matematicka funkce
 */
interface IFunction
{
    /**
	 * Vypocte a vrati hodnotu funkce v zadanem bode
	 */
    double ValueAt(double p);

    /// <summary>
    /// Vypočte a vrátí hodnotu derivace v bodě <paramref name="x"/>.
    /// </summary>
    /// <param name="x">Bod pro výpočet hodnoty derivace</param>
    /// <returns>Vrací hodnotu derivace v bodě <paramref name="x"/></returns>
    double Differentiate(double x);
}

/**
 * Trida pro numericky vypocet urciteho integralu funkce
 */
class Integrator
{
    /** Krok pro vypocet integralu */
    double delta;

    /**
	 * Numbericky vypocte a vrati urcity integral zadane fukce f od a do b
	 */
    public double Integrate(IFunction f, double a, double b)
    {
        double result = 0;
        double p = a;
        double v = f.ValueAt(p);    //výška y
        while (p + delta < b)
        {
            //obdelniky sirky delta
            result += delta * v;
            p += delta;
            v = f.ValueAt(p);
        }
        // jeste posledni obdelnik, ktery bude uzsi nez delta
        result += (b - p) * v;
        return result;
    }

    /**
	 * Nastavi krok pro vypocet integralu
	 * @param d krok pro vypocet integralu
	 */
    public void SetDelta(double d)
    {
        this.delta = d;
    }
}

/**
 * Linerani funkce
 * @author Libor Vasa
 */
class LinearFunction : AbstractFunction
{
    /** Smernice funkce */
    double k;
    /** Posun funkce */
    double q;

    /**
	 * Vytvori novou linearni funkci se zadanymi koeficienty
	 */
    public LinearFunction(double k, double q)
    {
        this.k = k;
        this.q = q;
    }

    /// <summary>
    /// Funkce vypočte funkční hodnotu lineární funkce v bodě <paramref name="p"/>. 
    /// </summary>
    /// <param name="p">Bod pro výpočet funkční hodnoty</param>
    /// <returns>Funkční hodnota v bodě <paramref name="p"/></returns>
    public override double ValueAt(double p)
    {
        return (k * p + q);
    }
}

/// <summary>
/// Třída reprezentující kvadratickou funkci
/// </summary>
class QuadraticPolynomial : AbstractFunction
{
    double a, b, c;
    /// <summary>
    /// Konstruktor třídy <b>QuadraticPolynomial</b> přebírající koeficienty funkce.
    /// </summary>
    /// <param name="a">Koeficient kvadratického členu</param>
    /// <param name="b">Koeficient lineárního členu</param>
    /// <param name="c">Koeficient absolutního členu</param>
    public QuadraticPolynomial(double a, double b, double c)
    {
        this.a = a;
        this.b = b;
        this.c = c;
    }

    /// <summary>
    /// Funkce vypočte funkční hodnotu kvadratické funkce v bodě <paramref name="x"/>.
    /// </summary>
    /// <param name="x">Bod pro výpočet funkční hodnoty</param>
    /// <returns>Funkční hodnota v bodě <paramref name="x"/></returns>
    public override double ValueAt(double x)
    {
        return (a * Math.Pow(x, 2) + b * x + c);
    }
}

/// <summary>
/// Abstraktní třída poskytující implementaci numerické derivace.
/// </summary>
abstract class AbstractFunction : IFunction
{
    double epsilon = 0.01;

    /// <summary>
    /// Funkce vrací hodnotu derivace funkce v bodě <paramref name="x"/>.
    /// </summary>
    /// <param name="x">Bod pro výpočet hodnoty derivace</param>
    /// <returns>Vrací hodnotu derivace funkce v bodě <paramref name="x"/></returns>
    public double Differentiate(double x)
    {
        double h = 0.1;
        double result1 = (ValueAt(x + h) - ValueAt(x))/h;
        h /= 2;
        double result2 = (ValueAt(x + h) - ValueAt(x)) / h;

        while(Math.Abs(result1-result2) > epsilon)
        {
            result1 = result2;
            h /= 2;
            result2 = (ValueAt(x + h) - ValueAt(x)) / h;
        }

        return result1;
    }

    /// <summary>
    /// Setter pro nastavení atributu epsilon.
    /// </summary>
    /// <param name="newEpsilon">Nová hodnota epsilon</param>
    public void SetEpsilon(double newEpsilon)
    {
        epsilon = newEpsilon;
    }

    /// <summary>
    /// Vypočte a vrátí hodnotu funkce v zadaném bodě.
    /// </summary>
    /// <param name="p">Bod pro výpočet hodnoty funkce</param>
    /// <returns>Vrací hodnotu funkce v bodě <paramref name="p"/></returns>
    public abstract double ValueAt(double p);
}

/// <summary>
/// Třída reprezentující obecný polynom.
/// </summary>
class GeneralPolynomial : AbstractFunction
{
    double[] coefs;
    /// <summary>
    /// Kontruktor třídy <b>GeneralPolynomial</b> přebírající pole reálných čísel představující koeficienty polynomu.
    /// </summary>
    /// <param name="coefs">Pole reálných čísel představující koeficienty polynomu</param>
    public GeneralPolynomial(double[] coefs)
    {
        this.coefs = coefs;
    }

    /// <summary>
    /// Funkce vrací hodnotu polynomu v bodě <paramref name="p"/>.
    /// </summary>
    /// <param name="p">Bod pro výpočet hodnoty polynomu</param>
    /// <returns>Vrací hodnotu polynomu v bodě <paramref name="p"/></returns>
    public override double ValueAt(double p)
    {
        double result = 0;
        for(int i = 0; i < coefs.Length; i++)
        {
            /*
            if (coefs.Length - i - 1 == 0)
            {
                result += coefs[i];
            }
            else
            {
                result += coefs[i] * Math.Pow(p, coefs.Length - i - 1);
            }
            */
            result += coefs[i] * Math.Pow(p, coefs.Length - i - 1);
        }
        return result;
    }
}

/// <summary>
/// Třída reprezentující funkci sinus.
/// </summary>
class Sine : AbstractFunction
{
    /// <summary>
    /// Funkce vypočítá hodnotu funkce sinus v bodě  <paramref name="p"/>.
    /// </summary>
    /// <param name="p">Bod pro výpočet hodnoty funkce sinus</param>
    /// <returns>Vrací hodnotu funkce v bodě <paramref name="p"/></returns>
    public override double ValueAt(double p)
    {
        return Math.Sin(p);
    }
}

/// <summary>
/// Třída reprezentující derivaci funkce pomocí numerické derivace.
/// </summary>
class Derivative : AbstractFunction
{
    IFunction f;

    /// <summary>
    /// Konstruktor třídy <b>Derivative</b> přebírající jako parametr funkci.
    /// </summary>
    /// <param name="f">Derivovaná funkce</param>
    public Derivative(IFunction f)
    {
        this.f = f;
    }

    /// <summary>
    /// Funkce vypočte hodnotu derivace v bodě <paramref name="p"/>.
    /// </summary>
    /// <param name="p">Bod pro výpočet hodnoty derivace</param>
    /// <returns>Vrací hodnotu derivace v bodě <paramref name="p"/></returns>
    public override double ValueAt(double p)
    {
        return f.Differentiate(p);
    }
}