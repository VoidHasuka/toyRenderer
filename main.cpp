#include "tgaimage.h"
#include "model.h"
#include "rasterizer.h"

#pragma region ������ʼ��

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const int width = 800;
const int height = 800;
const int depth = 255;

TGAImage image(width, height, TGAImage::RGB);
Model* model = new Model("obj/african_head.obj");
float* zbuffer = new float[width * height];

Vec3f cameraPos(0, 0, 3);//������ڷŵ�λ��
Vec3f lightDir(0, 0, -1);

#pragma endregion



//��
//void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
//	
//
//
//	//����һ
//	/*for (float t = 0.; t < 1.; t += .1) {
//		int x = x0 + (x1 - x0) * t;
//		int y = y0 + (y1 - y0) * t;
//		image.set(x, y, color);
//	}*/    //��һ���ֲܴڵ�ֱ��
//
//
//	//������
//	//for (int x = x0; x <= x1; x++) {
//	//	float t = (x - x0) / (float)(x1 - x0);  //��ǰ��Ⱦ�ٷֱ�
//	//	int y = y0 * (1. - t) + y1 * t;
//	//	image.set(x, y, color);
//	//}   //�����ܴ���������ǣ����x>x1���޷���Ⱦ������Ҳ����һ�����⣬�������̶ȹ̶�����ȡ����x����dx>dy��ʱ��Ƚ��ʺ�ʹ�ã���dx<dyʱ�Ͳ�������
//
//
//
//	//������
//	bool steep = false; //�Ƿ�Ϊ���͵��ߣ���dy>dxΪtrue����ʱ���Ǳ�֤��Ⱦ��˳���Ǵ����½ǿ�ʼ��t��ȡ����x������y��ֵ
//	if (abs(x0 - x1) < abs(y0 - y1))
//	{
//		steep = true;
//		std::swap(x0, y0);
//		std::swap(x1, y1);
//	}
//	if (x0 > x1) {
//		std::swap(x0, x1);
//		std::swap(y0, y1);
//	}
//	int dx = x1 - x0;  //�˴������Ż�������ѭ�������еĸ��������㣬ע���ʱ��ֵ����y�����
//	int dy = y1 - y0;
//	//float derror = abs(dy / (float)dx); //б��
//	//float error = 0;
//	//int y = y0;
//	//for (int x = x0; x <= x1; x++) {
//	//	if (steep) {
//	//		image.set(y, x, color);
//	//	}
//	//	else {
//	//		image.set(x, y, color);
//	//	}
//
//	//	error += derror; //�������л��ۣ������õ�y0+derror>y0+0.5��ȡy=y0+1����֮y=y0����y0+2*derror<y0+0.5����ȡy=y0
//	//	if (error > .5) {
//	//		y += (y1 > y0 ? 1 : -1);
//	//		error -= 1;
//	//	}
//	//}
//
//	/*
//	�Ż�ԭ���������б��Ϊdy/dx(slope)����x0��Ϊx0+1ʱ����Ϊ y0+1*slope < y0+0.5 ?y0:y0+1����ϸ�۲죬�����������Ż�Ϊ 2dy<dx ?y:y+1����ʱ����ȥ�˸������ļ���
//	xÿ�����ƶ�һ�Σ�error�ͻ���һ��2dy��ֵ��������error(���2dy)����dxʱ���� y += 1������ȥ2dx
//	Ϊʲô����ȥ2dx���˴�ȡ����һ������ĸ�����ǶԱȵ���dx��������Ǹպõ���dx��ô������һ��������0.5������ֱ�Ӽ���1������������Ȼ�͵÷�����
//	*/
//
//	int derror2 = abs(dy) * 2; //��һ���Ż���������ȥ����������
//	int error2 = 0;
//	int y = y0;
//	for (int x = x0; x <= x1; x++) {
//		if (steep) {
//			image.set(y, x, color);
//		}
//		else {
//			image.set(x, y, color);
//		}
//		error2 += derror2;
//		if (error2 > dx) {
//			y += (y1 > y0 ? 1 : -1);
//			error2 -= dx * 2;
//		}
//	}
//}

//��������ĳ�����������ڵ���������
//Vec3f barycentric(Vec3f* pts, Vec3f p) {
//	
//	float xa = pts[0].x;
//	float ya = pts[0].y;
//	float xb = pts[1].x;
//	float yb = pts[1].y;
//	float xc = pts[2].x;
//	float yc = pts[2].y;
//	float x = p.x;
//	float y = p.y;
//
//	float gamma =(float) ((ya - yb) * x + (xb - xa) * y + xa * yb - xb * ya) / ((ya - yb) * xc + (xb - xa) * yc + xa * yb - xb * ya);
//	float beta =(float) ((ya - yc) * x + (xc - xa) * y + xa * yc - xc * ya) / ((ya - yc) * xb + (xc - xa) * yb + xa * yc - xc * ya);
//	float alpha = 1 - gamma - beta;
//
//	return Vec3f(gamma, beta, alpha);
//}

//����������ת��Ϊƽ�����꣬ע�⣬�������������-1 ~ 1��������
//Vec3f worldToScrren(Vec3f v) {
//
//	return Vec3f(int((v.x + 1.) * width / 2. + .5), int((v.y + 1.) * height / 2. + .5), v.z);
//
//}

//������
//void triangle(Vec3f* pts,Vec2i*uvs,float* zBuffer, TGAImage& image, float intensity)
//{
//
//	/*line(t0.x, t0.y, t1.x, t1.y, image, green);
//	line(t1.x, t1.y, t2.x, t2.y, image, green);
//	line(t2.x, t2.y, t0.x, t0.y, image, red);*/
//
//
//	////�Զ����������ʹ��t1��t2��t3Ϊy��С����
//	//if (t0.y > t1.y)std::swap(t0, t1);
//	//if (t0.y > t2.y)std::swap(t0, t2);
//	//if (t1.y > t2.y)std::swap(t1, t2);
//
//	//int total_height = std::max( t2.y - t0.y,1);
//	//int middle_height =std::max(t1.y - t0.y,1);
//	//int small_height = std::max(t2.y - t1.y, 1);
//
//	////�����°벿��
//	//for (int t = t0.y; t <= t1.y; t++) {
//	//	
//	//	float alpha = (float)(t - t0.y) / (total_height);
//	//	float beta = (float)(t - t0.y) / (middle_height);
//
//	//	Vec2i A = t0 + (t2 - t0) * alpha;
//	//	Vec2i B = t0 + (t1 - t0) * beta;
//
//	//	image.set(A.x, t, color);
//	//	image.set(B.x, t, color);
//
//	//	if (A.x > B.x)std::swap(A, B);
//
//	//	//����ˮƽ��ɫ��
//	//	for (int x = A.x; x <= B.x; x++) {
//	//		image.set(x, t, color);
//	//	}
//	//}
//
//	////�����ϰ벿��
//	//for (int t = t1.y; t <= t2.y; t++) {
//	//	
//	//	float alpha = (float)(t - t0.y) / (total_height);
//	//	float beta = (float)(t - t1.y) / (small_height);
//
//	//	Vec2i A = t0 + (t2 - t0) * alpha;
//	//	Vec2i B = t1 + (t2 - t1) * beta;
//
//	//	image.set(A.x, t, color);
//	//	image.set(B.x, t, color);
//
//	//	if (A.x > B.x)std::swap(A, B);
//
//	//	//����ˮƽ��ɫ��
//	//	for (int x = A.x; x <= B.x; x++) {
//	//		image.set(x, t, color);
//	//	}
//
//	//}
//
//	//���л���
//	//for (int t = t0.y; t <= t2.y; t++) {
//	//	
//	//	bool half = t >= t1.y;
//
//	//	float alpha = (float)(t - t0.y) / (total_height);
//	//	//float beta = 0;
//	//	//if (t > t1.y) beta = (float)(t - t1.y) / (small_height);
//	//	//else  beta = (float)(t - t0.y) / (middle_height);
//	//	float beta = (float)(half ? t - t1.y : t - t0.y) / (half ? small_height : middle_height);
//
//
//	//	Vec2i A = t0 + (t2 - t0) * alpha;
//	//	//Vec2i B;
//	//	//if(t>t1.y)B = t1 + (t2 - t1) * beta;
//	//	//else B = t0 + (t1 - t0) * beta;;
//	//	Vec2i B = half ? (t1 + (t2 - t1) * beta) : (t0 + (t1 - t0) * beta);
//
//	//	image.set(A.x, t, color);
//	//	image.set(B.x, t, color);
//
//	//	if (A.x > B.x)std::swap(A, B);
//
//	//	//����ˮƽ��ɫ��
//	//	for (int x = A.x; x <= B.x; x++) {
//	//		image.set(x, t, color);
//	//	}
//	//}
//
//	/*
//	��Ϊʹ�ð�Χ������Ⱦ����ĳ�����ڸ��������ڲ�����Ⱦ
//	*/
//
//	//Vec2i bboxMin(image.get_width() - 1, image.get_height() - 1);
//	//Vec2i bboxMax(0, 0);
//
//	////���������εİ�Χ��
//	//bboxMin.x = std::min({ pts[0].x, pts[1].x, pts[2].x, bboxMin.x });
//	//bboxMin.y = std::min({ pts[0].y, pts[1].y, pts[2].y, bboxMin.y });
//	//bboxMax.x = std::max({ pts[0].x, pts[1].x, pts[2].x, bboxMax.x });
//	//bboxMax.y = std::max({ pts[0].y, pts[1].y, pts[2].y, bboxMax.y });
//
//	////������Χ�����������أ������������ڲ�����Ⱦ
//	//for (int x = bboxMin.x; x <= bboxMax.x; x++) {
//	//	for (int y = bboxMin.y; y <= bboxMax.y; y++) {
//	//		Vec2i P(x, y);
//	//		Vec3f barycentricP = barycentric(pts, P);
//
//	//		if (barycentricP.x < -0.01 || barycentricP.y < -0.01 || barycentricP.z < -0.01)
//	//			continue;
//
//	//		image.set(x, y, color);
//	//	}
//	//}
//
//	//����Z-Buffer��Ⱦ
//
//	float minx = std::min({ pts[0].x,pts[1].x,pts[2].x });
//	float maxx = std::max({ pts[0].x,pts[1].x,pts[2].x });
//	float miny = std::min({ pts[0].y,pts[1].y,pts[2].y });
//	float maxy = std::max({ pts[0].y,pts[1].y,pts[2].y });
//
//
//	int min_x = (int)std::floor(minx);  //����ȡ��
//	int max_x = (int)std::ceil(maxx);   //����ȡ��
//	int min_y = (int)std::floor(miny);
//	int max_y = (int)std::ceil(maxy);
//
//	//������Χ��������
//	for (int i = min_x; i <= max_x; i++) {
//		for (int j = min_y; j <= max_y; j++) {
//
//			Vec3f P(i, j, 0);
//			Vec2i uvP;
//			Vec3f barycentricP = barycentric(pts, P);
//
//			if (barycentricP.x < -0.01 || barycentricP.y < -0.01 || barycentricP.z < -0.01)continue;
//
//			//�����Ĳ�ֵ����������ε�zֵ  ����zBuffer���±��Ƕ�ά��һάת��������һ��(x,y)��Ӧ�±�Ϊ(x+width*y)
//			float z_interpolation = barycentricP.x * pts[0].z + barycentricP.y * pts[1].z + barycentricP.z * pts[2].z;
//			uvP = uvs[0] * barycentricP.x + uvs[1] * barycentricP.y + uvs[2] * barycentricP.z;
//			if (z_interpolation > zBuffer[static_cast<int>(i + j * width)]) {
//
//				zBuffer[static_cast<int>(i + j * width)] = z_interpolation;
//				TGAColor color = model->diffuse(uvP);
//				image.set(P.x, P.y, TGAColor(color.r*intensity,color.g*intensity,color.b*intensity,255));
//			}
//		}
//	}
//}


#pragma region ��ͼ�任
//
////���ֲ�����任Ϊ�������
//Matrix local2homo(Vec3f v) {
//	Matrix m(4, 1);
//	m[0][0] = v.x;
//	m[1][0] = v.y;
//	m[2][0] = v.z;
//	m[3][0] = 1.0f;
//
//	return m;
//}
//
////ģ�ͱ任���� δ��
//Matrix modelMatrix() {
//	return Matrix::identity(4);
//}
//
////��ͼ�任����  δ��
//Matrix viewMatrix() {
//	return Matrix::identity(4);
//}
//
////͸��ͶӰ�任����
//Matrix projectionMatrix() {
//	Matrix projection = Matrix::identity(4);
//	projection[3][2] = -1.0f / cameraPos.z;
//	return projection;
//}
//
////͸�ӳ���
//Matrix projectionDivision(Matrix m) {
//	m[0][0] = m[0][0] / m[3][0];
//	m[1][0] = m[1][0] / m[3][0];
//	m[2][0] = m[2][0] / m[3][0];
//	m[3][0] = 1.0f;
//	return m;
//}
//
////�ӿڱ任����NDC����ת��Ϊ��Ļ����
//Matrix viewportMatrix(int x, int y, int w, int h) {
//	Matrix m = Matrix::identity(4);
//	m[0][3] = x + w / 2.f;
//	m[1][3] = y + h / 2.f;
//	m[2][3] = depth / 2.f;
//
//	m[0][0] = w / 2.f;
//	m[1][1] = h / 2.f;
//	m[2][2] = depth / 2.f;
//	return m;
//}
//
////���������ָ�Ϊ��ά����
//Vec3f homo2vertices(Matrix m)
//{
//	return Vec3f(m[0][0], m[1][0], m[2][0]);
//}

#pragma endregion


//ģ�ͱ任����
Matrix modelMatrix()
{
    return Matrix::identity(4);
}

//��ͼ�任����
Matrix viewMatrix()
{
    return Matrix::identity(4);
}

//͸��ͶӰ�任����
Matrix projectionMatrix()
{
    Matrix projection = Matrix::identity(4);
    projection[3][2] = -1.0f / cameraPos.z;
    return projection;
}

Vec3f vertex_shader(const vertex_shader_payload& payload)
{
    return payload.position;
}

Vec3f normal_fragment_shader(const fragment_shader_payload& payload)
{
    Vec3f normal_frag = payload.normal;
    Vec3f return_color = (normal_frag.normalize() + Vec3f(1.0f, 1.0f, 1.0f)) * 0.5;

    return Vec3f(return_color.x * 255, return_color.y * 255, return_color.z * 255);
}

Vec3f F_fragment_shader(const fragment_shader_payload& payload)
{
    
    Vec3f color_frag = payload.color;
    return color_frag;
}   

Vec3f G_fragment_shader(const fragment_shader_payload& payload)
{
    Vec3f color_frag = payload.color;
    return color_frag;
}

int main(int argc, char** argv) {
    if (2 == argc) {
        model = new Model(argv[1]);
    }
    else {
        model = new Model("african_head.obj");
    }

    std::cout << model->nfaces() << "  " << model->nverts() << std::endl;

    //����TGAͼ��
    TGAImage image(width, height, TGAImage::Format::RGB);

    //�洢������Ҫ���Ƶ���������Ƭ

    //������դ������
    rst::rasterizer r(width, height);

    //������������
    Texture tex("african_head_diffuse.tga");
    r.set_texture(tex);

    //���֡�����zBuffer
    r.clear(rst::Buffers::Color);
    r.clear(rst::Buffers::Depth);

    //����MVP����
    r.set_model(modelMatrix());
    r.set_view(viewMatrix());
    r.set_projection(projectionMatrix());

    //���ö�����ɫ����ƬԪ��ɫ��
    r.set_vertexShader(vertex_shader);
    r.set_fragmentShader(normal_fragment_shader);

    //����ģ��
    r.draw(model->TriangleList);
    //��֡�����е���ɫֵд��image��
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            Vec3f color = r.frame_buffer[j * width + i];
            image.set(i, j, TGAColor(color.x, color.y, color.z, 255));
        }
    }
    image.flip_vertically();
    image.write_tga_file("output.tga");
}

