// Реализовать классы "полином", "полином Лагранжа", "полином Ньютона"

#include <iostream>
#include <iomanip>
#include <cstdlib>
using namespace std;
#include <math.h>
#define N 5 // Количество значений сетки
#define A -2 // Начало промежутка
#define B 1 // Конец промежутка

struct net; // Сетка значений
struct elem; // Элемент полинома
class polynom; // Полином
class polynom_lagrange; // Полином в форме Лагранжа
class polynom_newton; //  Полином в форме Ньютона

struct net // Сетка значений
{
	double x, y; 
	void init(net* n); // Функция, инициализирующая сетку значениями
	double count_y(double x) { return x*x*x+5*x; } // Вычисление у
	void print(net* n) const; // Печать сетки
};

void net::print(net* n) const // Печать сетки
{
	int i;
	for (i=0; i<N; i++)
		cout<<"x"<<i<<" = "<<setw(13)<< right <<n[i].x<<"	"<<"y"<<i<<" = "<< setw(13) << right << n[i].y<<endl;
}

void net::init(net* n) // Функция, инициализирующая сетку значениями
{
	double k=fabs(B-A)/N; // Шаг
	for (int i=0; i<N; i++)
	{
		if (i==0) // Записываем начало промежутка
			n[i].x=A;
		else
			if (i==N-1) // Записываем конец промежутка
				n[i].x=B;
			else // Или просто записываем значение, которое отличается от предыдущего на 1 шаг
				n[i].x=n[i-1].x+k; 
		n[i].y=count_y(n[i].x); // Вычисляем значение y
	}
}

struct elem // Элемент полинома
{
	elem(double k=1, int e=0, elem* n=0): exp(e), koef(k), next(n) {} // Конструктор
	double koef; // Коэффициент
	int exp; // Степень 
	elem* next; // Указатель на следующий элемент
	void print() const { cout<<koef<<"x^"<<exp<<" "; }
};

class polynom // Полином
{
public:
	friend class polynom_lagrange;
	friend class polynom_newton;
	elem* get_head() { return head; }
	void set_head(elem* el) { head=el; }
	polynom() { head=0; } // Конструктор по умолчанию
	polynom(double k, int e); // Конструктор, создающий моном для вызова из классов "полином Лагранжа" и "полином Ньютона"
	polynom(const polynom& p); // Копирующий конструктор
	~polynom() { del(head); } // Деструктор
	friend ostream& operator << (ostream& out, const polynom& p); // Вывод полинома
	polynom& operator = (const polynom& p); // Оператор присваивания
	polynom& operator += (const polynom& p); // Оператор += для полинома и полинома
	polynom& operator += (double num); // Оператор += для полинома и числа
	polynom operator + (const polynom& p); // Оператор сложения для двух полиномов
	polynom operator + (double num); // Оператор сложения полинома и числа
	friend polynom operator + (double num, const polynom& p); // Дружественная функция для сложения числа и полинома
	polynom operator * (const polynom& p); // Оператор умножения двух полиномов
	polynom operator * (double num); // Оператор умножения полинома на число
	polynom& operator *= (const polynom& p); // Оператор умножения с присваиванием для двух полиномов
	polynom& operator *= (double num); // Оператор умножения с присваиванием для полинома и числа
	polynom& operator -= (double num); // Оператор разности полинома и числа с присваиванием 
	friend polynom operator * (double num, const polynom& p); // Дружественная функция для умножения числа на полином
	elem* add_elem(elem* h, double k, int e); // Функция добавления элемента в конец списка
	double count_horner(elem* h, double x); // Вычисление значения функции в точке по правилу Горнера
	
private:
	elem* head;	// Голова списка
	elem* add(int b, elem* tmp1, elem* tmp2, elem* tmp3, elem* rez); // Добавление элемента при умножении
	elem* paste_elem(elem* el, double k, int e); // Вставка элемента между двумя другими
	void del(elem* h); // Удаление списка
	elem* del_elem(elem* h, elem* el); // Удаление элемента
	elem* copy(elem* t, elem* const h); // Функция копирования списка
	elem* sum(elem* const h1, elem* const h2); // Функция сложения
	elem* multiply(elem* const h1, elem* const h2); // Функция умножения
	void add_numeric(elem* h, double num); // Добавление числа
	void mult_numeric(elem* h, double num); // Умножение на число
	void pass_elem(double* y, double x, int count); // Функция для счёта в точке, пропускающая элементы с нулевыми коэфф-тами
};

polynom::polynom(double k, int e) // Конструктор, создающий моном для вызова из класса "полином Лагранжа"
{
	elem* el=new elem(k, e); 
	head=el;
} 

void polynom::mult_numeric(elem* h, double num) // Умножение на число
{
	elem* tmp=h;
	while (tmp!=0)
	{ // Пока не дошли до конца списка, умножаем каждый элемент на число
		tmp->koef*=num;
		tmp=tmp->next;
	}
}

void polynom::add_numeric(elem* h, double num) // Добавление числа
{
	elem* tmp=h; // Копия головы
	while (tmp->next!=0) // Идём до конца списка
		tmp=tmp->next;
	if (tmp->exp==0) // Если последний элемент - число, то просто прибавляем к нему новое число
		tmp->koef+=num;
	else // Иначе добавляем число в конец списка как элемент с нулевой степенью
		h=add_elem(h, num, 0);	
}

double polynom::count_horner(elem* h, double x) // Вычисление значения функции в точке по правилу Горнера
{
	elem* temp=h; // Копия головы
	double y; // Будущий результат вычисления
	if (temp->next==0) // Если в полиноме только 1 элемент:
	{  
		y=temp->koef*x; // Умножаем 1ый раз
		pass_elem(&y, x, temp->exp); // Домножаем оставшееся кол-во раз
		return y;
	}	
	y=h->koef; // Заносим в результат коэффициент при старшей степени
	while(temp->next!=0)
	{ // Пока не закончится список из элементов полинома, считаем значение в точке х
		y*=x;
		pass_elem(&y, x, temp->exp - temp->next->exp); // Пропускаем элементы с нулевыми коэффициентами, домножая на х
		y+=temp->next->koef; // Прибавляем коэффициент при следующем элементе
		temp=temp->next;	
	}
	pass_elem(&y, x, temp->exp);	
	if (temp->exp!=0)
		y*=x; // Если последним стоит элемент с ненулевой степенью, нужно домножить последний раз на х
	return y;
}

void polynom::pass_elem(double* y, double x, int count) // Функция для счёта в точке, пропускающая элементы с нулевыми коэфф-тами
{
	while(count>1)
	{
		(*y)*=x;
		count--;
	}
}

elem* polynom::paste_elem(elem* el, double k, int e) // Вставка элемента между двумя другими
{
	elem* n=new elem(k, e); // Создаём элемент с нужными коэффициентами
	elem* tmp=el; // Копия элемента, после которого нужно вставить новый
	elem* tmp1=el->next; // Элемент, который должен стоять после нового
	tmp->next=n;
	n->next=tmp1; 
	return tmp1;
}

elem* polynom::del_elem(elem* h, elem* el) // Удаление элемента
{
	elem* temp=h; 
	while(temp->next!=el)
		temp=temp->next; // Идём по списку, пока не найдём элемент, который нужно удалить
	temp->next=el->next; // Предыдущий элемент должен ссылаться на следующий после удаляемого
	delete el; // Удаляем сам элемент
	return temp;
}

elem* polynom::add(int b, elem* tmp1, elem* tmp2, elem* tmp3, elem* rez) // Добавление элемента при умножении
{
	switch(b)
	{
		case 0:
		{ // Добавление элемента в конец:
			tmp3=add_elem(rez, tmp1->koef*tmp2->koef, tmp1->exp+tmp2->exp);
			break;	
		}
		case 1:
		{ // При равных степенях меняем коэффициент при элементе или удаляем его, если сумма коэффициентов даёт 0:
			if (((tmp1->koef*tmp2->koef) + tmp3->koef) == 0)
				return del_elem(rez, tmp3); 
			else
				tmp3->koef+=tmp1->koef*tmp2->koef;
			break;
		}
		case 2: // Вставка элемента м/у 2мя другими:
			tmp3=paste_elem(tmp3, (tmp1->koef*tmp2->koef), (tmp1->exp+tmp2->exp));
	}
	return tmp3;
}

elem* polynom::add_elem(elem* h, double k, int e) // Функция добавления элемента в конец списка
{
	elem* tmp=h;
	elem* n=new elem(k, e); // Создается элемент с нужными параметрами
	if (h==0) // Для пустого списка назначаем голову
		h=n; 
	else // Если в списке уже есть элементы, то присоединяем последний
	{
		while (tmp->next!=0)
			tmp=tmp->next;
		tmp->next=n;
	}
	return h;
}

polynom::polynom(const polynom& p) // Копирующий конструктор
{
	head=0;
	head=copy(head, p.head);
}

void polynom::del(elem* h) // Удаление списка
{
	if (h!=0){ // Если список не пуст, происходит его удаление
		elem* temp=h;
		while (temp!=0 )
		{ // Пока не дойдём до конца списка	
			elem* n=temp->next; // Запоминаем следующий элемент
			delete temp; // Удаляем текущий элемент
			temp=n; // Переходим в следующему элементу
		}
		h=0; // Обнуляем указатель на голову списка
	}
}

elem* polynom::multiply(elem* const h1, elem* const h2) // Функция умножения двух списков
{
	elem* tmp1=h1; // Начало 1го списка
	elem* tmp2=h2; // Начало 2го списка
	elem* rez=0; // Новый список
	elem* tmp3=0; // Начало нового списка
	int count=0, b=0;
	while (tmp1!=0)
	{
		while (tmp2!=0) 
		{
			while (count!=0 && tmp3!=0) 
			{ // Проходим по новому списку, если уже умножили хотя бы 1 раз
				if (tmp1->exp+tmp2->exp==tmp3->exp) // Если у будущего элемента степень совпадает со степенью эл-та из нового списка
				{
					b=1;
					break;
				}
				else
				{
					if (tmp3->next!=0) // Если в новом списке ещё есть элементы
					{
						if (tmp1->exp+tmp2->exp < tmp3->exp && tmp1->exp+tmp2->exp > tmp3->next->exp) 
						{
							b=2; // Случай, когда элемент можно вставить м/у 2мя другими
							break;
						}
					}
				}
				tmp3=tmp3->next; 
			}
			tmp3=add(b, tmp1, tmp2, tmp3, rez); // Добавляем новый элемент
			tmp2=tmp2->next; 
			if (rez==0) rez=tmp3; // Если список ещё пуст, назначаем ему голову
			b=0;
		}
		tmp1=tmp1->next;
		count++; // Считаем, сколько раз умножили
		tmp2=h2; // Встаём на начало 2го списка
		tmp3=rez; // Встаём на начало нового списка
	}
	return rez; // Возвращается голова сформированного списка
}

elem* polynom::sum(elem* const h1, elem* const h2) // Функция сложения двух списков
{
	elem* tmp1=h1; // Начало 1го списка
	elem* tmp2=h2; // Начало 2го списка
	elem* n=0; // Новый список
	elem* h=0; // Голова нового списка
	do
	{ // Пока не закончится один из списков, добавляем в 1ый список элементы в нужной последовательности
		if (tmp1->exp==tmp2->exp) // Если у элементов одна и та же степень:
		{
			n=add_elem(n, (tmp1->koef+tmp2->koef), tmp1->exp); // Добавляем эл-т, полученный в рез-тате сложения
			tmp1=tmp1->next; // Переходим к следующему элементу 1го списка
			tmp2=tmp2->next; // Переходим к следующему элементу 2го списка
		}
		else 
		{ // Добавляется тот элемент, который имеет большую степень
			if (tmp1->exp>tmp2->exp)
			{ 
				n=add_elem(n, tmp1->koef, tmp1->exp);
				tmp1=tmp1->next;
			}
			else
			{
				n=add_elem(n, tmp2->koef, tmp2->exp);
				tmp2=tmp2->next;
			}
		}
		if (h==0) h=n; // Если был добавлен 1ый элемент, нужно назначить голову списка
	}
	while (tmp1!=0 && tmp2!=0);	
	if (tmp2!=0) // Если кончился 1ый список, то остаток 2го копируется в новый	
		n=copy(n, tmp2); 
	if (tmp1!=0) // Если кончился 2ой список, то остаток 1го копируется в новый
		n=copy(n, tmp1); 
	return h;
}

elem* polynom::copy(elem* t, elem* const h) 
{ // Функция копирования списка (t-конец 1го списка, h-элемент 2го списка, с которого начинается копирование)
	if (h==0) return t;
	elem* h1=h; // Начало 2го списка
	elem* t1=0; // Начало присоединяемой части в 1ом списке
	if (h!=0)
	{
		do
		{
			t1=add_elem(t1, h1->koef, h1->exp);
			h1=h1->next; 
			if (t1==0) t1=t; // Если был добавлен 1ый элемент, нужно назначить голову списка
		} while (h1!=0); // Пока не дойдём до последнего элемента	
	}
	return t1;
}

/* Операторы */

ostream& operator << (ostream& out, const polynom& p) // Вывод элементов списка
{
	if (p.head==0)
		out<<0<<endl; // Если полином пуст, выводим 0
	else
	{
		elem* temp=p.head;
		do
		{ // Печатаем элементы полинома, пока не дойдём до конца
			if (temp->koef>0 && temp!=p.head)
				out<<"+ "; // Если коэффициент положительный и при этом не первый в списке, то перед печатью ставим +
			temp->print();
			temp=temp->next;
		}
		while (temp!=0); 
	}
	out<<endl;
	return out;
}

polynom& polynom::operator = (const polynom& p) // Оператор присваивания
{
	del(head); // Очищаем полином
	head=copy(head, p.head); // Копируем полином
	return *this;
}

polynom& polynom::operator += (const polynom& p) // Оператор += для полинома и полинома
{ 
	polynom n;
	if (head==0)
	{ // Если полином, к которому прибавляем, равен нулю, то просто копируем в него второй полином
		head=copy(head, p.head); 
		return *this;
	}
	if (p.head==0) // Если второй полином равен нулю, то ничего не прибавляем
		return *this; 
	n.head=sum(head, p.head);
	del(head);
	head=copy(head, n.head);
	return *this;
}

polynom& polynom::operator += (double num) // Оператор += для полинома и числа
{
	if (head==0) // Если полином пуст, то сразу добавляем новый элемент
		head=add_elem(head, num, 0);
	else
		add_numeric(head, num);
	return *this;
}

polynom polynom::operator + (const polynom& p) // Оператор сложения для двух полиномов
{
	polynom n;
	if (head!=0 && p.head!=0) // Если оба полинома не равны 0, то складываем их
		n.head=sum(head, p.head);
	else
	{ // Если один из полиномов равен 0, то копируем другой в новый полином
		if (head==0) 
			n.head=copy(n.head, p.head); 
		else
			n.head=copy(n.head, head);
	}
	return n;
}

polynom polynom::operator + (double num) // Оператор сложения полинома и числа
{
	polynom n;
	n.head=copy(n.head, head);
	if (n.head==0) // Если полином пуст, то сразу добавляем новый элемент
		n.head=add_elem(n.head, num, 0);
	else
		add_numeric(n.head, num);
	return n;
}

polynom operator + (double num, const polynom& p) // Дружественная функция для сложения числа и полинома
{
	polynom n(p);
	return n+=num;
}

polynom polynom::operator * (const polynom& p) // Оператор умножения двух полиномов
{
	polynom n;
	if (head==0 || p.head==0) // Если хотя бы один из полиномов равен нулю, возвращаем нулевой полином
		return n; 
	n.head=multiply(head, p.head);
	return n;
} 

polynom polynom::operator * (double num) // Оператор умножения полинома на число
{
	polynom n(*this);
	mult_numeric(n.head, num);
	return n;
}

polynom& polynom::operator *= (double num) // Оператор умножения с присваиванием для полинома и числа
{
	if (head==0) // Если полином пуст, сразу возвращаем его 
		return *this;
	if (num==0)
	{ // Если число - ноль, то очищаем полином и возвращаем пустоту
		del(head);
		return *this;
	}
	mult_numeric(head, num);
	return *this;
}

polynom& polynom::operator *= (const polynom& p) // Оператор умножения с присваиванием для двух полиномов
{
	if (head==0) // Если один из полиномов равен 0, то возвращаем нулевой полином
		return *this;
	if (p.head==0)
	{
		del(head);
		return *this;
	}
	polynom n;
	n.head=multiply(head, p.head); // Заносим результат умножения в новый полином
	del(head); // Удаляем старый полином
	head=copy(head, n.head); // Копируем в него только что созданный полином
	return *this;
} 

polynom& polynom::operator -= (double num) // Оператор разности с присваиванием 
{
	return *this += (-num);
}

polynom operator * (double num, const polynom& p) // Дружественная функция для умножения числа на полином
{
	polynom n;
	if (num==0 || p.head==0) // Если число равно нулю или полином пустой, возвращаем 0 в кач-ве результата
		return n;
	n.head=n.copy(n.head, p.head); // Иначе копируем в новый полином старый и возвращаем результат умножения на число
	return n*=num;
}

class polynom_lagrange // Полином в форме Лагранжа
{
public:
	polynom_lagrange(net* n) { pol=lagrange(n); } // Конструктор
	friend ostream& operator << (ostream& out, const polynom_lagrange& l) { out<<l.pol; return out; }
	polynom& get_pol() { return pol; }
private:
	polynom pol;
	polynom lagrange(net* n); // Функция, формирующая полином в форме Лагранжа
};

polynom polynom_lagrange::lagrange(net* n) // Функция, формирующая полином в форме Лагранжа
{
	polynom l;
	for (int i=0; i<N; i++)
	{
		polynom p; // Числитель 
		double z=1; // Знаменатель
		for (int j=0; j<N; j++)
		{
			polynom bin(1, 1); // Создаём моном 
			if (i!=j)
			{
				bin-=n[j].x; // Вычитаем число, чтобы получить бином
				if (p.head==0)
					p=bin; // Если числитель ещё пуст, записываем в него созданный бином
				else p*=bin; // Иначе домножаем на очередной бином
				z*=(n[i].x - n[j].x); // Считаем знаменатель, домножая на очередную разность иксов
			}
		}
		l+=(n[i].y / z) * p; // К результату прибавляем полученное слагаемое
	}
	return l;
}

class polynom_newton // Полином в форме Ньютона
{
public:
	polynom_newton(net* n) { pol=newton(n); } // Конструктор
	friend ostream& operator << (ostream& out, const polynom_newton& n) { out<<n.pol; return out; }
	polynom& get_pol() { return pol; }
private:
	polynom pol;
	polynom newton(net* n); // Функция, формирующая полином в форме Ньютона
	void diff(int i, net* n, polynom& r); // Формирование разделённой разности
	void mult(int i, net* n, polynom& r); // Функция для домножения разделённых разностей на биномы (x-x(0)) ... (x-x(i-1))
};

void polynom_newton::diff(int i, net* n, polynom& r) // Формирование разделённой разности
{
	for (int j=0; j<=i; j++)
	{
		double z=1; // Знаменатель разделённой разности
		for (int k=0; k<=i; k++)
		{
			if (k!=j) // Умножаем разности иксов, формируя знаменатель
				z*=n[j].x-n[k].x; 
		}
		r+=n[j].y / z; // Добавляем к результату очередное слагаемое, формируя разделённую разность
	}
}

void polynom_newton::mult(int i, net* n, polynom& r) 
{ // Функция для домножения разделённых разностей на биномы (x-x(0)) ... (x-x(i-1))
	for (int k=0; k<i; k++)
	{ 
		polynom bin(1, 1);
		bin-=n[k].x;
		r*=bin;
	}
}

polynom polynom_newton::newton(net* n) // Функция, формирующая полином в форме Ньютона
{
	polynom nt;
	nt+=n[0].y; // Сначала полином равен разделённой разности нулевого порядка
	for (int i=1; i<N; i++)
	{
		polynom r; 
		diff(i, n, r); // Формирование разделённой разности
		mult(i, n, r); // Домножение разделённых разностей на биномы (x-x(0)) ... (x-x(i-1))
		nt+=r; 
	}
	return nt;
}

double count(polynom& p, double n) // Функция счёта в точке (общая)
{
	return p.count_horner(p.get_head(), n);
}

void print_pair(polynom_lagrange& l, polynom_newton& nt, net* n, double x)
{ // Функция для вывода значений х и у
	cout << "x=" << setw(14) << right << x << "	y(regular)=" << setw(14) << right << n->count_y(x)
	<<"	y(Lagrange)=" << setw(14) << right << count(l.get_pol(), x) << "	y(Newton)="
	<< setw(14) <<right << count(nt.get_pol(), x) << endl;
}

void print(polynom_lagrange& l, polynom_newton& nt, net* n) // Вывод таблицы значений 
{
	double half=(n[1].x-n[0].x)/2; // Половина шага
	double x;
	for (int i=0; i<N; i++)
	{
		x=n[i].x; // Выводим значения в узлах сетки
		print_pair(l, nt, n, x);
		// Выводим промежуточные значения
		x+=half;
		print_pair(l, nt, n, x);
	}
}

int main()
{
	polynom p, r, q, d;
	cout.setf(ios::scientific);
	p.set_head(p.add_elem(p.get_head(), 7, 4));
	p.set_head(p.add_elem(p.get_head(), -14, 2));
	p.set_head(p.add_elem(p.get_head(), 8, 1));
	r.set_head(r.add_elem(r.get_head(), 2, 3));
	r.set_head(r.add_elem(r.get_head(), 7, 1));
	cout<<"Given data set:"<<endl;
	net* n=new net[N];
	n->init(n);
	n->print(n);
	polynom_lagrange l(n);
	cout<<"\nLagrange polynomial:\n"<<l;
	polynom_newton nt(n);
	cout<<"\nNewton polynomial:\n"<<nt<<endl;
	cout<<"Values in and between the given points of interpolation:"<<endl;
	print(l, nt, n);
	delete [] n;
	system("pause");
}
