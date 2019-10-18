using namespace std;

#define GRAY unsigned char

struct COLOR{
	unsigned char r,g,b;
};

template <class datatype> class Image
{
	public:
		
		int W;			//幅
		int H;			//高さ
		datatype **data;//画像配列

		
		void load( char *filename );	 //画像ロード
		void create( int x, int y );	 //配列確保
		void copy( Image<datatype> img );//配列コピー
		void save( char *filename );	 //画像保存
		
		Image(){data=NULL;};
		~Image(){
			if(data != NULL){
				for( int i = 0; i < H; i++)
					delete[] data[i];

				delete[] data;
			}
		}

		Image( const Image &img ){
			W = img.W;
			H = img.H;
			data = allocate();
			for(int i = 0; i < H; i++ )
				memcpy( data[i], img.data[i], sizeof(datatype) * W );
		}

	private:
		datatype **allocate()		//メモリ確保
		{
			datatype **d = new datatype*[H];

			for( int i=0; i < H; i++ )
				d[i] = new datatype[W];

			return d;
		}
};


/******************************************************************************/
/***							画像読み込み								***/
/******************************************************************************/
template <class datatype> void Image<datatype>::load( char *filename )
{
	FILE *fp;
	char str[256];
	
	if(( fp = fopen( filename, "rb" ) ) == NULL ){
  	  puts("File not open!ファイルが存在しない可能性があります.");
		exit(1);
	}
	else
		printf("FileLoad: %s \n", filename );

	fgets( str, 100, fp ); 	//P5,P6スキャン
	if( strcmp( str, "P6\n" ) != 0 &&  strcmp( str, "P5\n" ) != 0 ){	//読み込む画像がppmかpgmでなければ強制終了
		puts(" THIS FILE CANNOT USE! ファイルがppm,pgmではありません");
		exit(1);
	}

	//エラー処理(読み込む画像と配列の型があっているか判断)
	if( sizeof(datatype)==sizeof(COLOR) && strcmp( str, "P6\n" ) != 0 ){
				puts("COLOR型はpgmファイルを読み込めません．");
				exit(1);
	}
	else if( sizeof(datatype)==sizeof(GRAY) && strcmp( str, "P5\n" ) != 0 ){
				cout<<"GRAY型はppmファイルを読み込めません．"<<endl;
				exit(1);
	}
	
	do{									//コメント行読み飛ばし
		fgets( str, 100, fp ); 	
	}while( strncmp( str, "#", 1 ) == 0 );

	sscanf( str, "%d %d\n", &W, &H );	//画像サイズ読み込み
	fgets( str, 100, fp );				//階調数読みとばし

	create( W, H ); //data確保
	
	for( int i=0; i < H; i++ ){			//ロード
		for( int j=0; j < W; j++ ){
			fread( &data[i][j], sizeof(datatype), 1, fp );
		}
	}

	fclose(fp);
}

/******************************************************************************/
/***						空の画像を生成する								***/
/******************************************************************************/
template <class datatype> void Image<datatype>::create( int x, int y )
{
	W = x;
	H = y;
	data = allocate();
}

/******************************************************************************/
/***						画像のコピー									***/
/******************************************************************************/
template <class datatype> void Image<datatype>::copy( Image<datatype> img )
{
	for( int i = 0; i < img.H; i++ ){
		for( int j = 0; j < img.W; j++ ){
			data[i][j] = img.data[i][j];
		}
	}
}

/******************************************************************************/
/***								保存									***/
/******************************************************************************/
template <class datatype> void Image<datatype>::save( char *filename )
{
	FILE *fp = fopen( filename, "wb" );
	
	//ヘッダ書き込み
	if( sizeof(datatype) == sizeof(COLOR) ){
		fprintf( fp, "P6\n%d %d\n255\n", W, H); 
	}
	else if( sizeof(datatype) == sizeof(GRAY) ){
			fprintf( fp, "P5\n%d %d\n255\n", W, H);
	}
	else{
		printf("この型は保存できません\n");
	}
	//画像保存
	for( int i = 0; i < H; i++ ){			
		for( int j = 0; j < W; j++ ){
			fwrite( &data[i][j], sizeof(datatype), 1, fp );
		}
	}
	
	fclose(fp);

	printf("File Save : %s \n", filename );	//セーブファイル名表示
}
