ImageIO使用方法

GRAY	: グレースケールの画像
COLOR	: カラーの画像

●Imageの中身
	W 	: 画像の幅  (pixel)
	H 	: 画像の高さ(pixel)
	data 	: 画像の二次元配列 data[0][0] ～ data[H-1][W-1]
	
		COLOR型の場合　	data[j][i].r	<- (i,j)番目の画素の赤色の輝度値
				data[j][i].g	<- (i,j)番目の画素の緑色の輝度値
				data[j][i].b	<- (i,j)番目の画素の青色の輝度値

		GRAY型の場合　	data[j][i]	<- (i,j)番目の画素の輝度値


①配列の宣言：

	書き方: Image<datatype> name;

		datatype : 画像配列の型を指定 ( GRAYかCOLOR )
		name     : 画像配列の名前　(好きな名前でOK)

	(例) Image<COLOR> img;  
		 imgという名前のCOLOR型の配列ができる．

②画像の読み込み：
	
	書き方：name.load( image name );

		load 	   : "画像を読み込め！"という命令文
		image name : 読み込みたい画像の名前
		
	(例) img.load( "lenna.ppm" );
		lenna.ppmという名前の画像をimg配列に読み込む
	
    (※注)
	・画像を読み込むと,W(画像の幅),H(画像の高さ)は自動的に画像サイズに応じた
	　値が代入される．
	　((例)320x240の画像を読み込んだ場合，Wに320,Hに240が代入される．)	

	・読み込める画像は".ppm"または".pgm"拡張子のみ．".bmp"等は
	　読み込めないので注意．
	
	・宣言した配列の型と画像の型が異なる場合読み込めない．
	　たとえば配列の宣言を"GRAY"型とし，".ppm"のファイル(COLOR型)を
	  読み込もうとしたときエラーが発生する．
　　　　　(例)
		Image<GRAY> gray;
		gray.load( "lenna.ppm" );  <- 配列と画像の型が違うのでエラー

③画像の保存：
	
	書き方: name.save( save name );
		
		save	   : "画像を保存せよ！"という命令文
		save name  : 保存する画像の名前
	
	(例) img.save( "output.ppm" );
		img配列を output.ppm という名前をつけて保存


④配列のみの宣言（画像を読み込まず，配列だけほしい場合）

	書き方：Image<datatype> name;　<- まず宣言
		name.create( width, height );
			
		create	: "画像配列を作りなさい"という命令文
		width	: 画像の幅(pixel)	 
		height	: 画像の高さ(pixel)

	(例) Image<COLOR> output;
	     output.create( 320, 240 );
		320x240のサイズの画像配列をCOLOR型で作成．

⑤画像配列のコピー
	
	書き方：name.copy( img );
		
		copy	: "画像配列imgをnameにコピーしなさい"という命令
		img	: コピー元の画像配列

	(例) Image<COLOR> img,output;
	　　 img.load( "lenna.ppm" );		<- lenna.ppmをロード
	     output.create( img.W, img.H );  	<- lenna.ppm　と同じ画像サイズの
						    配列作成
	     output.copy ( img ); 		<- outputにimgがコピーされる

　　(※注)
	・違う型同士のコピーはできない．
	  
	　(例)	Image<COLOR> img;
		Image<GRAY> output;	　	
		img.load( "lenna.ppm" );	<- lenna.ppmをロード
	     	output.create( img.W, img.H );  <- lenna.ppm　と同じ画像サイズの
						    配列作成
		output.copy ( img ); 		<- コピーする配列が違う型なのでエラー！！
