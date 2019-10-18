using namespace std;

#define GRAY unsigned char

struct COLOR{
	unsigned char r,g,b;
};

template <class datatype> class Image
{
	public:
		
		int W;			//��
		int H;			//����
		datatype **data;//�摜�z��

		
		void load( char *filename );	 //�摜���[�h
		void create( int x, int y );	 //�z��m��
		void copy( Image<datatype> img );//�z��R�s�[
		void save( char *filename );	 //�摜�ۑ�
		
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
		datatype **allocate()		//�������m��
		{
			datatype **d = new datatype*[H];

			for( int i=0; i < H; i++ )
				d[i] = new datatype[W];

			return d;
		}
};


/******************************************************************************/
/***							�摜�ǂݍ���								***/
/******************************************************************************/
template <class datatype> void Image<datatype>::load( char *filename )
{
	FILE *fp;
	char str[256];
	
	if(( fp = fopen( filename, "rb" ) ) == NULL ){
  	  puts("File not open!�t�@�C�������݂��Ȃ��\��������܂�.");
		exit(1);
	}
	else
		printf("FileLoad: %s \n", filename );

	fgets( str, 100, fp ); 	//P5,P6�X�L����
	if( strcmp( str, "P6\n" ) != 0 &&  strcmp( str, "P5\n" ) != 0 ){	//�ǂݍ��މ摜��ppm��pgm�łȂ���΋����I��
		puts(" THIS FILE CANNOT USE! �t�@�C����ppm,pgm�ł͂���܂���");
		exit(1);
	}

	//�G���[����(�ǂݍ��މ摜�Ɣz��̌^�������Ă��邩���f)
	if( sizeof(datatype)==sizeof(COLOR) && strcmp( str, "P6\n" ) != 0 ){
				puts("COLOR�^��pgm�t�@�C����ǂݍ��߂܂���D");
				exit(1);
	}
	else if( sizeof(datatype)==sizeof(GRAY) && strcmp( str, "P5\n" ) != 0 ){
				cout<<"GRAY�^��ppm�t�@�C����ǂݍ��߂܂���D"<<endl;
				exit(1);
	}
	
	do{									//�R�����g�s�ǂݔ�΂�
		fgets( str, 100, fp ); 	
	}while( strncmp( str, "#", 1 ) == 0 );

	sscanf( str, "%d %d\n", &W, &H );	//�摜�T�C�Y�ǂݍ���
	fgets( str, 100, fp );				//�K�����ǂ݂Ƃ΂�

	create( W, H ); //data�m��
	
	for( int i=0; i < H; i++ ){			//���[�h
		for( int j=0; j < W; j++ ){
			fread( &data[i][j], sizeof(datatype), 1, fp );
		}
	}

	fclose(fp);
}

/******************************************************************************/
/***						��̉摜�𐶐�����								***/
/******************************************************************************/
template <class datatype> void Image<datatype>::create( int x, int y )
{
	W = x;
	H = y;
	data = allocate();
}

/******************************************************************************/
/***						�摜�̃R�s�[									***/
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
/***								�ۑ�									***/
/******************************************************************************/
template <class datatype> void Image<datatype>::save( char *filename )
{
	FILE *fp = fopen( filename, "wb" );
	
	//�w�b�_��������
	if( sizeof(datatype) == sizeof(COLOR) ){
		fprintf( fp, "P6\n%d %d\n255\n", W, H); 
	}
	else if( sizeof(datatype) == sizeof(GRAY) ){
			fprintf( fp, "P5\n%d %d\n255\n", W, H);
	}
	else{
		printf("���̌^�͕ۑ��ł��܂���\n");
	}
	//�摜�ۑ�
	for( int i = 0; i < H; i++ ){			
		for( int j = 0; j < W; j++ ){
			fwrite( &data[i][j], sizeof(datatype), 1, fp );
		}
	}
	
	fclose(fp);

	printf("File Save : %s \n", filename );	//�Z�[�u�t�@�C�����\��
}
