ImageIO�g�p���@

GRAY	: �O���[�X�P�[���̉摜
COLOR	: �J���[�̉摜

��Image�̒��g
	W 	: �摜�̕�  (pixel)
	H 	: �摜�̍���(pixel)
	data 	: �摜�̓񎟌��z�� data[0][0] �` data[H-1][W-1]
	
		COLOR�^�̏ꍇ�@	data[j][i].r	<- (i,j)�Ԗڂ̉�f�̐ԐF�̋P�x�l
				data[j][i].g	<- (i,j)�Ԗڂ̉�f�̗ΐF�̋P�x�l
				data[j][i].b	<- (i,j)�Ԗڂ̉�f�̐F�̋P�x�l

		GRAY�^�̏ꍇ�@	data[j][i]	<- (i,j)�Ԗڂ̉�f�̋P�x�l


�@�z��̐錾�F

	������: Image<datatype> name;

		datatype : �摜�z��̌^���w�� ( GRAY��COLOR )
		name     : �摜�z��̖��O�@(�D���Ȗ��O��OK)

	(��) Image<COLOR> img;  
		 img�Ƃ������O��COLOR�^�̔z�񂪂ł���D

�A�摜�̓ǂݍ��݁F
	
	�������Fname.load( image name );

		load 	   : "�摜��ǂݍ��߁I"�Ƃ������ߕ�
		image name : �ǂݍ��݂����摜�̖��O
		
	(��) img.load( "lenna.ppm" );
		lenna.ppm�Ƃ������O�̉摜��img�z��ɓǂݍ���
	
    (����)
	�E�摜��ǂݍ��ނ�,W(�摜�̕�),H(�摜�̍���)�͎����I�ɉ摜�T�C�Y�ɉ�����
	�@�l����������D
	�@((��)320x240�̉摜��ǂݍ��񂾏ꍇ�CW��320,H��240����������D)	

	�E�ǂݍ��߂�摜��".ppm"�܂���".pgm"�g���q�̂݁D".bmp"����
	�@�ǂݍ��߂Ȃ��̂Œ��ӁD
	
	�E�錾�����z��̌^�Ɖ摜�̌^���قȂ�ꍇ�ǂݍ��߂Ȃ��D
	�@���Ƃ��Δz��̐錾��"GRAY"�^�Ƃ��C".ppm"�̃t�@�C��(COLOR�^)��
	  �ǂݍ������Ƃ����Ƃ��G���[����������D
�@�@�@�@�@(��)
		Image<GRAY> gray;
		gray.load( "lenna.ppm" );  <- �z��Ɖ摜�̌^���Ⴄ�̂ŃG���[

�B�摜�̕ۑ��F
	
	������: name.save( save name );
		
		save	   : "�摜��ۑ�����I"�Ƃ������ߕ�
		save name  : �ۑ�����摜�̖��O
	
	(��) img.save( "output.ppm" );
		img�z��� output.ppm �Ƃ������O�����ĕۑ�


�C�z��݂̂̐錾�i�摜��ǂݍ��܂��C�z�񂾂��ق����ꍇ�j

	�������FImage<datatype> name;�@<- �܂��錾
		name.create( width, height );
			
		create	: "�摜�z������Ȃ���"�Ƃ������ߕ�
		width	: �摜�̕�(pixel)	 
		height	: �摜�̍���(pixel)

	(��) Image<COLOR> output;
	     output.create( 320, 240 );
		320x240�̃T�C�Y�̉摜�z���COLOR�^�ō쐬�D

�D�摜�z��̃R�s�[
	
	�������Fname.copy( img );
		
		copy	: "�摜�z��img��name�ɃR�s�[���Ȃ���"�Ƃ�������
		img	: �R�s�[���̉摜�z��

	(��) Image<COLOR> img,output;
	�@�@ img.load( "lenna.ppm" );		<- lenna.ppm�����[�h
	     output.create( img.W, img.H );  	<- lenna.ppm�@�Ɠ����摜�T�C�Y��
						    �z��쐬
	     output.copy ( img ); 		<- output��img���R�s�[�����

�@�@(����)
	�E�Ⴄ�^���m�̃R�s�[�͂ł��Ȃ��D
	  
	�@(��)	Image<COLOR> img;
		Image<GRAY> output;	�@	
		img.load( "lenna.ppm" );	<- lenna.ppm�����[�h
	     	output.create( img.W, img.H );  <- lenna.ppm�@�Ɠ����摜�T�C�Y��
						    �z��쐬
		output.copy ( img ); 		<- �R�s�[����z�񂪈Ⴄ�^�Ȃ̂ŃG���[�I�I
