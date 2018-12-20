#include<bits/stdc++.h>
using namespace std;

const int N=62536;
typedef bitset<N+1> poly;
poly pw[N*2],f;
poly mul(const poly &a,const poly &b){
	poly c;c.reset();
	for(int i=0;i<N;i++)if(a[i]){
		for(int j=0;j<N;j++)if(b[j]){
			c^=pw[i+j];
		}
	}
	return c;
}
poly two[N*2];

int deg(const poly &p){
	int n=0;
	for(int i=N;i>=0;i--)if(p[i]){
		n=i;
		return i;
	}
	return -1e9;
}

namespace FFT{
	#define MAXN 262144
	#define DB double
	const DB pi = acos(-1);
	struct CP {
		DB x, y; CP(){} inline CP(DB a, DB b):x(a),y(b){}
		inline CP operator + (const CP&r) const { return CP(x + r.x, y + r.y); }
		inline CP operator - (const CP&r) const { return CP(x - r.x, y - r.y); }
		inline CP operator * (const CP&r) const { return CP(x*r.x-y*r.y, x*r.y+y*r.x); }
		inline CP conj() { return CP(x, -y); }
	} a[MAXN], b[MAXN], t;
	int n, m;
	void FFT(CP*a, int n, int f) {
		int i, j, k;
		for(i = j = 0; i < n; ++ i) {
			if(i > j) swap(a[i], a[j]);
			for(k = n>>1; (j ^= k) < k; k >>= 1);
		}
		for(i = 1; i < n; i <<= 1) {
			CP wn(cos(pi/i), f*sin(pi/i));
			for(j = 0; j < n; j += i<<1) {
				CP w(1, 0);
				for(k = 0; k < i; ++ k, w=w*wn) {
					CP x = a[j+k], y = w*a[i+j+k];
					a[j+k] = x+y; a[i+j+k] = x-y;
				}
			}
		}
		if(-1 == f) for(i = 0; i < n; ++ i) a[i].x /= n;
	}
	poly mul(const poly &p1,const poly &p2){
		poly ans;
		n = deg(p1); m = deg(p2);
		if(n<0||m<0)return ans;
		for(int i = 0; i <= n; ++ i) a[i].x = p1[i];
		for(int i = 0; i <= m; ++ i) a[i].y = p2[i];
		for(m += n, n = 1; n <= m; n <<= 1);
		FFT(a, n, 1); CP qua(0, -0.25);
		for(int i = 0, j; i < n; ++ i) j = (n-i)&(n-1), b[i] = (a[i]*a[i]-(a[j]*a[j]).conj())*qua;
		FFT(b, n,-1);
		for(int i = 0; i <= m; ++ i){
			int x=int(b[i].x+0.2);
			if(x%2)ans^=pw[i];
		}
		for(int i=0;i<n;i++)a[i]=b[i]=CP();
		return ans;
	}
	
}



unsigned char XTIME(unsigned char x) {  
    return ((x << 1) ^ ((x & 0x80) ? 0x1b : 0x00));  
}  
unsigned char mul(unsigned char a, unsigned char b) {  
    unsigned char temp[8] = { a };  
    unsigned char tempmultiply = 0x00;  
    int i = 0;  
    for (i = 1; i < 8; i++) {  
        temp[i] = XTIME(temp[i - 1]);  
    }  
    tempmultiply = (b & 0x01) * a;  
    for (i = 1; i <= 7; i++) {  
        tempmultiply ^= (((b >> i) & 0x01) * temp[i]);  
    }  
    return tempmultiply;  
}  


poly mod(poly a,poly m){
	int n=deg(m);	
	poly r;
	for(int i=N;i>=0;i--){
		r=r<<1;
		if(a[i])
			r[0]=1;
		if(r[n])
			r^=m;
	}
	return r;
}
poly zero,one;
poly gcd(poly a,poly b){
	int da=deg(a);
	int db=deg(b);
	if(da>db){
		swap(da,db);
		swap(a,b);
	}
	int mx=db;
	double st=clock();
	while(1){

		if(a==zero)return b;
		int da=deg(a);
		int db=deg(b);
		double has=clock()-st;
		double all=has*mx/(mx-da+1);
		printf("deg: %5d %5d %.4fs\r",da,db,all-has);
		
		
		b=mod(b,a);
		swap(a,b);
	}
	puts("");
}
int last=1;
double st;
void get(int n){
	
	printf("get %d\n",n);
	two[0][1]=1;
	for(int i=last;i<=n;i++){
		double all=2*N*(clock()-st)/CLOCKS_PER_SEC/i;
		printf("%5d / %5d  %.4fs\r",i,N*2,all-(clock()-st)/CLOCKS_PER_SEC);
		two[i]=FFT::mul(two[i-1],two[i-1]);
	}
	puts("");
	last=n+1;
}

bool is_irr(poly f){
	int n=deg(f);
	poly g,h;
	vector<int>prime;
	int t=n;
	for(int i=2;i*i<=n;i++){
		if(n%i==0){
			prime.push_back(i);
			while(n%i==0)n/=i;
		}
	}
	if(n>1)prime.push_back(n);
	n=t;

	//cerr<<"first check"<<endl;
	//cerr<<"check size"<<prime.size()<<endl;
	reverse(prime.begin(),prime.end());
	last=1;
	for(auto p:prime){
		//cerr<<"check "<<p<<endl;
		int ni=n/p;
		get(ni);
		h=two[ni]^two[0];
		g=gcd(f,h);
	//	cerr<<"fgh"<<endl;
	//	cerr<<f.to_string()<<endl;	
	//	cerr<<g.to_string()<<endl;	
	//	cerr<<h.to_string()<<endl;	
		if(g==one){
		
		}else{
			return false;
		}
	}
	get(n);
	g=two[n]^two[0];
	if(g==zero){
	
	}else
		return false;
	return true;
}


void test_mul(){
	int t=1e6;
	while(t--){
		int x=rand()&255;
		int y=rand()&255;

		int res1=(int)mul(x,y);
		poly a,b;
		a=x;b=y;
		poly c=mul(a,b);
		assert(c.to_ulong()==res1);
	}
	cerr<<"test mul pass"<<endl;
}
void test_gcd(){
	int t=1e5;
	poly a=9;
	poly b=0x3f;
	poly g=gcd(a,b);
	cout<<g.to_ulong()<<endl;
	cerr<<"test mod pass"<<endl;
}

void setf(poly _f){
	f=_f;
	int n=deg(f);
	one[0]=1;
	pw[0][0]=1;
	for(int i=1;i<N*2;i++){
		pw[i]=pw[i-1]<<1;
		if(pw[i][n]){
			pw[i]^=f;
		}
	}


}
void test_irr(){
	for(int i=1;i<64;i++){
		poly f=i;
		setf(f);
		if(is_irr(f))
			cout<<f.to_string()<<endl;
	}
}
bool isprime(int n){
	for(int i=2;i*i<=n;i++)if(n%i==0)return false;
	return true;
}
int main(){
	st=clock();
	int start=N-(time(0)*351%2000);
	//cin>>start; 
	for(int i=0;i<100;i++){
		for(int j=1;j<10;j++){

			double st=clock();
			poly f;
			f[start-i]=1;
			f[j]=1;
			f[0]=1;
			
			int n=start-i;
			int k=j;
			//Swan's theorem
			
			if(n%2==0&&k%2==1){
				if(n*k/2%4==0 || n*k/2%4==1)
					continue;
			}
			if(n%2==1&&k%2==0){
				if(n%(2*k)&&(n%8==3||n%8==5))
					continue;
			}
			if(n%2==1&&k%2==0){
				if(n%(2*k)==0&&(n%8==1||n%8==7))
					continue;
			}
			
			if(!isprime(n))continue;
			cerr<<"test "<<start-i<<" "<<j<<endl;
			setf(f);
			cerr<<"after preprocess"<<endl;
			cout<<(clock()-st)/CLOCKS_PER_SEC<<endl;
			int t=is_irr(f);
			if(t){
				freopen("OwO.txt","w",stdout);
				cout<<"pass "<<start-i<<" "<<j<<endl;
				fclose(stdout);
				exit(0);
			}
		}

	}
	return 0;
}
