#ifdef GL_ES
precision mediump float;
#endif

uniform float time;
uniform vec2 mouse;
uniform vec2 resolution;

const float far = 10000.0;
vec2 intersectionRaySphereR(in vec3 ro, in vec3 rd, in vec3 center, const float radius2) {
    vec3 vdif = ro - center;
    float dt = dot(rd, vdif);
    float x = dt*dt - dot(vdif, vdif) + radius2;
    vec2 ret;
    if (x >= 0.0) {
        float sqt = sqrt(x);
        ret = vec2(-dt - sqt, -dt + sqt);
    }
    else {
        ret = vec2(far);
    }
    return ret;
}

vec3 mod289(vec3 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 mod289(vec4 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 permute(vec4 x) {
     return mod289(((x*34.0)+1.0)*x);
}

vec4 taylorInvSqrt(vec4 r)
{
  return 1.79284291400159 - 0.85373472095314 * r;
}

float snoise(vec3 v, out vec3 gradient)
{
  const vec2  C = vec2(1.0/6.0, 1.0/3.0) ;
  const vec4  D = vec4(0.0, 0.5, 1.0, 2.0);

// First corner
  vec3 i  = floor(v + dot(v, C.yyy) );
  vec3 x0 =   v - i + dot(i, C.xxx) ;

// Other corners
  vec3 g = step(x0.yzx, x0.xyz);
  vec3 l = 1.0 - g;
  vec3 i1 = min( g.xyz, l.zxy );
  vec3 i2 = max( g.xyz, l.zxy );

  //   x0 = x0 - 0.0 + 0.0 * C.xxx;
  //   x1 = x0 - i1  + 1.0 * C.xxx;
  //   x2 = x0 - i2  + 2.0 * C.xxx;
  //   x3 = x0 - 1.0 + 3.0 * C.xxx;
  vec3 x1 = x0 - i1 + C.xxx;
  vec3 x2 = x0 - i2 + C.yyy; // 2.0*C.x = 1/3 = C.y
  vec3 x3 = x0 - D.yyy;      // -1.0+3.0*C.x = -0.5 = -D.y

// Permutations
  i = mod289(i); 
  vec4 p = permute( permute( permute( 
             i.z + vec4(0.0, i1.z, i2.z, 1.0 ))
           + i.y + vec4(0.0, i1.y, i2.y, 1.0 )) 
           + i.x + vec4(0.0, i1.x, i2.x, 1.0 ));

// Gradients: 7x7 points over a square, mapped onto an octahedron.
// The ring size 17*17 = 289 is close to a multiple of 49 (49*6 = 294)
  float n_ = 0.142857142857; // 1.0/7.0
  vec3  ns = n_ * D.wyz - D.xzx;

  vec4 j = p - 49.0 * floor(p * ns.z * ns.z);  //  mod(p,7*7)

  vec4 x_ = floor(j * ns.z);
  vec4 y_ = floor(j - 7.0 * x_ );    // mod(j,N)

  vec4 x = x_ *ns.x + ns.yyyy;
  vec4 y = y_ *ns.x + ns.yyyy;
  vec4 h = 1.0 - abs(x) - abs(y);

  vec4 b0 = vec4( x.xy, y.xy );
  vec4 b1 = vec4( x.zw, y.zw );

  //vec4 s0 = vec4(lessThan(b0,0.0))*2.0 - 1.0;
  //vec4 s1 = vec4(lessThan(b1,0.0))*2.0 - 1.0;
  vec4 s0 = floor(b0)*2.0 + 1.0;
  vec4 s1 = floor(b1)*2.0 + 1.0;
  vec4 sh = -step(h, vec4(0.0));

  vec4 a0 = b0.xzyw + s0.xzyw*sh.xxyy ;
  vec4 a1 = b1.xzyw + s1.xzyw*sh.zzww ;

  vec3 p0 = vec3(a0.xy,h.x);
  vec3 p1 = vec3(a0.zw,h.y);
  vec3 p2 = vec3(a1.xy,h.z);
  vec3 p3 = vec3(a1.zw,h.w);

//Normalise gradients
  vec4 norm = taylorInvSqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
  p0 *= norm.x;
  p1 *= norm.y;
  p2 *= norm.z;
  p3 *= norm.w;

// Mix final noise value
  vec4 m = max(0.5 - vec4(dot(x0,x0), dot(x1,x1), dot(x2,x2), dot(x3,x3)), 0.0);
  vec4 m2 = m * m;
  vec4 m4 = m2 * m2;
  vec4 pdotx = vec4(dot(p0,x0), dot(p1,x1), dot(p2,x2), dot(p3,x3));

// Determine noise gradient
  vec4 temp = m2 * m * pdotx;
  gradient = -8.0 * (temp.x * x0 + temp.y * x1 + temp.z * x2 + temp.w * x3);
  gradient += m4.x * p0 + m4.y * p1 + m4.z * p2 + m4.w * p3;
  gradient *= 105.0;

  return 105.0 * dot(m4, pdotx);
}


void main( void ) {
	vec2 position = gl_FragCoord.xy;
	vec2 aspect = vec2(1.,resolution.y/resolution.x );
	vec2 p = 2.0*(position - 0.5*resolution) / resolution;
	vec2 m = vec2(2.0, 1.0)*(mouse - vec2(0.5,0.0))*vec2(0.25*3.14159, 0.5*3.14159)+vec2(0.0,-3.14159/16);
	vec3 u = vec3(0.0,1.0,0.0);
	vec3 v = vec3(aspect.y,0.0,0.0);
	vec3 ro = vec3(0.0, 0.0, -10.0);
	vec3 x = vec3(0.0,0.0,-9.0) + ro + p.x * u + p.y * v;
	vec3 rd = normalize(x - ro);
	vec3 cc = vec3(0.0, 0.0, 2.0);
	mat3 my = mat3(cos(m.y), 0.0, sin(m.y),
			0.0, 1.0, 0.0,
			-sin(m.y), 0.0, cos(m.y)
			);	 
	mat3 mx = mat3(1.0, 0.0, 0.0,
			0.0, cos(m.x), -sin(m.x),
			0.0, sin(m.x), cos(m.x)
			);	 

	mat3 myi = mat3(cos(-m.y), 0.0, sin(-m.y),
			0.0, 1.0, 0.0,
			-sin(-m.y), 0.0, cos(-m.y)
			);	 
	mat3 mxi = mat3(1.0, 0.0, 0.0,
			0.0, cos(-m.x), -sin(-m.x),
			0.0, sin(-m.x), cos(-m.x)
			);	 


	vec2 t = intersectionRaySphereR(ro, rd, cc, 0.45);
	vec3 ip = ro + t.x * rd;


	vec3 cc2 = cc + (mx*my*vec3(0.0,0.0,sqrt(0.45)));

	vec3 n = normalize(ip - cc);
	vec3 lpos = n+vec3(1.0,0.0,0.0);//+vec3(mouse,0.0);
	vec2 t2 = intersectionRaySphereR(ro, rd, cc2, 0.012);
	vec3 ip2 = ro + t2.x * rd;
	vec3 n2 = normalize(ip2 - cc2);
	vec3 col = vec3(0.5,0.2,0.2);
	vec3 g;
	float noise0 = sin(snoise(mxi*myi*11.0*vec3(ip-cc), g));
	bool isDvorec = distance(ip, cc2) < 0.25-0.015*noise0;
	vec3 dvorec = isDvorec ? vec3(0.7, 0.65, 0.65) : vec3(1.0);
	if(t2.x < t.x) {
		dvorec = vec3(0.8,0.75,0.75);
		t.x = t2.x;
		n = normalize(mix(n2, vec3(0.0, 0.0,1.0), 0.25*(1.0+noise0)));
		col = 0.75*vec3(0.65,0.3,0.35);
		ip = ip2;
	}
	float mult0 = isDvorec ? 18.0 : 1.0;
	float mult1 = isDvorec ? 1.0 : 3.0;
	float strength = isDvorec ? 0.23 : 0.1;
	float noise = 1.0 - strength + strength*(1.0+sin(mult1*snoise(mxi*myi*mult0*vec3(ip-cc), g)));

	float diff = t.x < far ? dot(normalize(ip - lpos), n) : 0.0;
	gl_FragColor=vec4(col*diff*dvorec*noise, 1.0);
}
