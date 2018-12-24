#include<bits/stdc++.h>
using namespace std;

const int N=8200;
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
void setf(poly _f){
	f=_f;
	int n=deg(f);
	pw[0][0]=1;
	for(int i=1;i<N*2;i++){
		pw[i]=pw[i-1]<<1;
		if(pw[i][n]){
			pw[i]^=f;
		}
	}
}
double mu=0.5;

double fac(int n){
	double ans=1;
	for(int i=2;i<=n;i++)ans*=n;
	return ans;
}
double C(int n,int m){
	return fac(n)/fac(m)/fac(n-m);
}

vector<double> fun(){
	poly g;
	g[7]=g[1]=g[0]=1;
	setf(g);
	double E=0;
	vector<double>ans;
	ans.resize(7);
	for(int i=0;i<(1<<7);i++){
		for(int j=0;j<(1<<7);j++){
			poly a=i,b=j;
			poly c=mul(a,b);
			double p=pow(mu,a.count())*pow(1-mu,7-a.count())*pow(mu,b.count())*pow(1-mu,7-b.count());
			for(int k=0;k<7;k++)if(c[k])
				ans[k]+=p;
		}
	}
	return ans;
}

double prob_cir[N*2];
double prob_line[N*2];

double bf_cir(int n){
	double ans=0;
	for(int S=0;S<(1<<n);S++){
		vector<int>a;
		for(int i=0;i<n;i++)
			a.push_back(S>>i&1);
		int res=0;
		double p=1;
		for(int i=0;i<n;i++){
			res^=a[i]&a[(i+1)%n];
			if(a[i])p*=mu;else p*=1-mu;
		}
		if(res==0)
			ans+=p;
	}
	return ans;
}
double bf_line(int n){
	double ans=0;
	for(int S=0;S<(1<<n);S++){
		vector<int>a;
		for(int i=0;i<n;i++)
			a.push_back(S>>i&1);
		int res=0;
		double p=1;
		for(int i=0;i+1<n;i++){
			res^=a[i]&a[(i+1)%n];
		}
		for(int i=0;i<n;i++){
			if(a[i])p*=mu;else p*=1-mu;
		}
		if(res==0)
			ans+=p;
	}
	return ans;
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
	void norm(vector<double>&p){
		double s=accumulate(p.begin(),p.end(),0.0);
		for(auto &x:p)x/=s;
	}
	vector<double> mul(const vector<double> &p1,const vector<double> &p2){
		vector<double> ans;
		n = p1.size(); m = p2.size();
		if(n<0||m<0)return ans;
		for(int i = 0; i <= n; ++ i) a[i].x = p1[i];
		for(int i = 0; i <= m; ++ i) a[i].y = p2[i];
		for(m += n, n = 1; n <= m; n <<= 1);
		ans.resize(m);
		FFT(a, n, 1); CP qua(0, -0.25);
		for(int i = 0, j; i < n; ++ i) j = (n-i)&(n-1), b[i] = (a[i]*a[i]-(a[j]*a[j]).conj())*qua;
		FFT(b, n,-1);
		for(int i = 0; i <= m; ++ i){
			ans[i]=b[i].x;
		}
		for(int i=0;i<n;i++)a[i]=b[i]=CP();
		return ans;
	}
	
}


double prob(int n,int m,int k){
	//\sum c_ix^i mod x^n+x^m+1
	// Pr[ x_k 's cof = 1 ]
	// c_k + c_{k+n-m} + c_{k+n}
	//cerr<<"prob "<<n<<" "<<m<<" "<<k<<endl;
	static int vis[N*2],deg[N*2];
	static vector<int> G[N*2];
	memset(vis,0,sizeof vis);
	for(int i=0;i<N*2;i++)G[i].clear();
	
	for(int i=0;i<n;i++){
		if(k-i>=0 && k-i<n){
			G[i].push_back(k-i+N);
			G[k-i+N].push_back(i);
		}
		if(k+n-m-i>=0 && k+n-m-i<n){
			G[i].push_back(k+n-m-i+N);
			G[k+n-m-i+N].push_back(i);
		}
		if(k+n-i>=0 && k+n-i<n){
			G[i].push_back(k+n-i+N);
			G[k+n-i+N].push_back(i);
		}
	}
	vector<double>pro;
	for(int i=0;i<n;i++){
		if(vis[i])continue;
		int u=i,len=0;
		bool is_cir=true;
		while(1){
			if(vis[u])break;
//			if(u<N)cerr<<"a"<<u;else cerr<<"b"<<u-N;cerr<<" ";
			vis[u]=1;len++;
			is_cir&=(G[u].size()==2);
			for(auto v:G[u]){
				if(!vis[v]){
					u=v;
				}
			}
		}
//		cerr<<endl;
//		cerr<<is_cir<<" "<<len<<endl;
		if(is_cir)
			pro.push_back(prob_cir[len]	);
		else
			pro.push_back(prob_line[len]);
	} 
	
	
	double dp[2]={1,0};
	double ans=0;
	for(auto p: pro){
		double x=p*dp[0]+(1-p)*dp[1];
		double y=(1-p)*dp[0]+p*dp[1];
		dp[0]=x;
		dp[1]=y;
	}
	
	return dp[1];
}

int main(){
	mu=1.0/405;
	
	static double dp[N*2][2][2];// 1..i determined a[i] is j ,Pr[Xor]=k;
	
	//a[1]=0
	
	double pro[2]={1-mu,mu};
	
	memset(dp,0,sizeof dp);
	dp[2][0][0]=1-mu;
	dp[2][1][0]=mu;
	for(int i=3;i<N*2;i++){
		for(int j=0;j<2;j++)
		for(int k=0;k<2;k++)
		for(int l=0;l<2;l++)
			dp[i][l][k^(j&l)]+=pro[l]*dp[i-1][j][k];
	}
	
	for(int i=0;i<N*2;i++){
		prob_cir[i]+=(1-mu)*(dp[i][0][0]+dp[i][1][0]);
		prob_line[i]+=(1-mu)*(dp[i][0][0]+dp[i][1][0]);
	}
	//a[1]=1
	
	memset(dp,0,sizeof dp);
	
	dp[2][0][0]=1-mu;
	dp[2][1][1]=mu;
	for(int i=3;i<N*2;i++){
		for(int j=0;j<2;j++)
		for(int k=0;k<2;k++)
		for(int l=0;l<2;l++)
			dp[i][l][k^(j&l)]+=pro[l]*dp[i-1][j][k];
	}
	
	for(int i=0;i<N*2;i++){
		prob_cir[i]+=mu*(dp[i][0][0]+dp[i][1][1]);
		prob_line[i]+=mu*(dp[i][0][0]+dp[i][1][0]);
	}
	
	//auto res=fun();
	//for(auto x:res)cout<<x<<",";cout<<endl;
	
	static vector<vector<double> > poly;
	int n=8100,m=9;
	poly.resize(n);
	for(int i=0;i<n;i++){
		double p=prob(n,m,i);
		cerr<<p<<endl;
		poly[i]=vector<double>{1-p,p};
	}
	while(poly.size()>1){
		int m=poly.size()/2;
		for(int i=0;i<m;i++){
			poly[i]=FFT::mul(poly[i],poly.back());
			FFT::norm(poly[i]);
			poly.pop_back();
		}
	}
	double res=0;
	for(int i=0;i<n;i+=10){
		for(int j=i;j<i+10&&j<n;j++)
			res+=poly[0][j];
		cout<<i<<" "<<res<<endl;
	}
	cout<<endl;

		
	/*for(int i=0;i<2*7-1;i++){
		printf("c_{%d}x^{%d}",i,i);
		if(i+1==2*7-1)printf("mod g=\n");
		else printf("+");
	}
	for(int i=0;i<7;i++){
		vector<int>cof;
		for(int j=0;j<2*7-1;j++)if(pw[j][i])
			cof.push_back(j);
		printf("(c_{%d}",cof[0]);
		for(int j=1;j<cof.size();j++)
			printf("+c_{%d}",cof[j]);
		printf(")x^{%d}",i);
		
		if(i+1==7)printf("\n");
		else printf("+");
	}*/
	
	
	
	return 0;
}
