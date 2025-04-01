#include "tgaimage.h"
#include "model.h"
#include "rasterizer.h"

#pragma region 变量初始化

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const int width = 800;
const int height = 800;
const int depth = 255;

TGAImage image(width, height, TGAImage::RGB);
Model* model = new Model("obj/african_head.obj");
float* zbuffer = new float[width * height];

Vec3f cameraPos(0, 0, 3);//摄像机摆放的位置
Vec3f lightDir(0, 0, -1);

#pragma endregion



//线
//void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
//	
//
//
//	//方法一
//	/*for (float t = 0.; t < 1.; t += .1) {
//		int x = x0 + (x1 - x0) * t;
//		int y = y0 + (y1 - y0) * t;
//		image.set(x, y, color);
//	}*/    //是一个很粗糙的直线
//
//
//	//方法二
//	//for (int x = x0; x <= x1; x++) {
//	//	float t = (x - x0) / (float)(x1 - x0);  //当前渲染百分比
//	//	int y = y0 * (1. - t) + y1 * t;
//	//	image.set(x, y, color);
//	//}   //并不能处理的问题是，如果x>x1将无法渲染，并且也存在一个问题，即采样程度固定，仅取决于x，当dx>dy的时候比较适合使用，但dx<dy时就不好用了
//
//
//
//	//方法三
//	bool steep = false; //是否为陡峭的线，即dy>dx为true，此时我们保证渲染的顺序是从左下角开始，t仅取决于x，仅对y插值
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
//	int dx = x1 - x0;  //此处进行优化，减少循环过程中的浮点数计算，注意此时插值仍在y上完成
//	int dy = y1 - y0;
//	//float derror = abs(dy / (float)dx); //斜率
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
//	//	error += derror; //对误差进行积累，如计算得到y0+derror>y0+0.5则取y=y0+1，反之y=y0，若y0+2*derror<y0+0.5则仍取y=y0
//	//	if (error > .5) {
//	//		y += (y1 > y0 ? 1 : -1);
//	//		error -= 1;
//	//	}
//	//}
//
//	/*
//	优化原理解析：记斜率为dy/dx(slope)，当x0进为x0+1时，则为 y0+1*slope < y0+0.5 ?y0:y0+1，仔细观察，该条件可以优化为 2dy<dx ?y:y+1，此时便消去了浮点数的计算
//	x每向右移动一次，error就积累一次2dy的值，当满足error(多个2dy)大于dx时则令 y += 1，并消去2dx
//	为什么是消去2dx？此处取用了一个镜像的概念，我们对比的是dx，如果我们刚好等于dx那么就是上一个方法的0.5，但我直接加了1，翻倍加了自然就得翻倍减
//	*/
//
//	int derror2 = abs(dy) * 2; //进一步优化，彻底消去浮点数计算
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

//三角形与某点在三角形内的重心坐标
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

//将世界坐标转换为平面坐标，注意，世界坐标局限于-1 ~ 1的正方体
//Vec3f worldToScrren(Vec3f v) {
//
//	return Vec3f(int((v.x + 1.) * width / 2. + .5), int((v.y + 1.) * height / 2. + .5), v.z);
//
//}

//三角形
//void triangle(Vec3f* pts,Vec2i*uvs,float* zBuffer, TGAImage& image, float intensity)
//{
//
//	/*line(t0.x, t0.y, t1.x, t1.y, image, green);
//	line(t1.x, t1.y, t2.x, t2.y, image, green);
//	line(t2.x, t2.y, t0.x, t0.y, image, red);*/
//
//
//	////对顶点进行排序，使得t1，t2，t3为y从小到大
//	//if (t0.y > t1.y)std::swap(t0, t1);
//	//if (t0.y > t2.y)std::swap(t0, t2);
//	//if (t1.y > t2.y)std::swap(t1, t2);
//
//	//int total_height = std::max( t2.y - t0.y,1);
//	//int middle_height =std::max(t1.y - t0.y,1);
//	//int small_height = std::max(t2.y - t1.y, 1);
//
//	////绘制下半部分
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
//	//	//绘制水平颜色线
//	//	for (int x = A.x; x <= B.x; x++) {
//	//		image.set(x, t, color);
//	//	}
//	//}
//
//	////绘制上半部分
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
//	//	//绘制水平颜色线
//	//	for (int x = A.x; x <= B.x; x++) {
//	//		image.set(x, t, color);
//	//	}
//
//	//}
//
//	//进行汇总
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
//	//	//绘制水平颜色线
//	//	for (int x = A.x; x <= B.x; x++) {
//	//		image.set(x, t, color);
//	//	}
//	//}
//
//	/*
//	改为使用包围盒来渲染，即某点若在该三角形内部则渲染
//	*/
//
//	//Vec2i bboxMin(image.get_width() - 1, image.get_height() - 1);
//	//Vec2i bboxMax(0, 0);
//
//	////计算三角形的包围盒
//	//bboxMin.x = std::min({ pts[0].x, pts[1].x, pts[2].x, bboxMin.x });
//	//bboxMin.y = std::min({ pts[0].y, pts[1].y, pts[2].y, bboxMin.y });
//	//bboxMax.x = std::max({ pts[0].x, pts[1].x, pts[2].x, bboxMax.x });
//	//bboxMax.y = std::max({ pts[0].y, pts[1].y, pts[2].y, bboxMax.y });
//
//	////遍历包围盒内所有像素，若在三角形内部则渲染
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
//	//采用Z-Buffer渲染
//
//	float minx = std::min({ pts[0].x,pts[1].x,pts[2].x });
//	float maxx = std::max({ pts[0].x,pts[1].x,pts[2].x });
//	float miny = std::min({ pts[0].y,pts[1].y,pts[2].y });
//	float maxy = std::max({ pts[0].y,pts[1].y,pts[2].y });
//
//
//	int min_x = (int)std::floor(minx);  //向下取整
//	int max_x = (int)std::ceil(maxx);   //向上取整
//	int min_y = (int)std::floor(miny);
//	int max_y = (int)std::ceil(maxy);
//
//	//遍历包围盒内像素
//	for (int i = min_x; i <= max_x; i++) {
//		for (int j = min_y; j <= max_y; j++) {
//
//			Vec3f P(i, j, 0);
//			Vec2i uvP;
//			Vec3f barycentricP = barycentric(pts, P);
//
//			if (barycentricP.x < -0.01 || barycentricP.y < -0.01 || barycentricP.z < -0.01)continue;
//
//			//用重心插值计算该三角形的z值  其中zBuffer的下标是二维向一维转化，即任一点(x,y)对应下标为(x+width*y)
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


#pragma region 视图变换
//
////将局部坐标变换为齐次坐标
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
////模型变换矩阵 未做
//Matrix modelMatrix() {
//	return Matrix::identity(4);
//}
//
////视图变换矩阵  未做
//Matrix viewMatrix() {
//	return Matrix::identity(4);
//}
//
////透视投影变换矩阵
//Matrix projectionMatrix() {
//	Matrix projection = Matrix::identity(4);
//	projection[3][2] = -1.0f / cameraPos.z;
//	return projection;
//}
//
////透视除法
//Matrix projectionDivision(Matrix m) {
//	m[0][0] = m[0][0] / m[3][0];
//	m[1][0] = m[1][0] / m[3][0];
//	m[2][0] = m[2][0] / m[3][0];
//	m[3][0] = 1.0f;
//	return m;
//}
//
////视口变换，将NDC坐标转换为屏幕坐标
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
////将齐次坐标恢复为三维坐标
//Vec3f homo2vertices(Matrix m)
//{
//	return Vec3f(m[0][0], m[1][0], m[2][0]);
//}

#pragma endregion


//模型变换矩阵
Matrix modelMatrix()
{
    return Matrix::identity(4);
}

//视图变换矩阵
Matrix viewMatrix()
{
    return Matrix::identity(4);
}

//透视投影变换矩阵
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

    //创建TGA图像
    TGAImage image(width, height, TGAImage::Format::RGB);

    //存储所有需要绘制的三角形面片

    //创建光栅化对象
    rst::rasterizer r(width, height);

    //给定纹理并设置
    Texture tex("african_head_diffuse.tga");
    r.set_texture(tex);

    //清空帧缓冲和zBuffer
    r.clear(rst::Buffers::Color);
    r.clear(rst::Buffers::Depth);

    //设置MVP矩阵
    r.set_model(modelMatrix());
    r.set_view(viewMatrix());
    r.set_projection(projectionMatrix());

    //设置顶点着色器和片元着色器
    r.set_vertexShader(vertex_shader);
    r.set_fragmentShader(normal_fragment_shader);

    //绘制模型
    r.draw(model->TriangleList);
    //将帧缓冲中的颜色值写入image中
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

